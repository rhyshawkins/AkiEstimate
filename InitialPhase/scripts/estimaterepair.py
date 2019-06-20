
import numpy
import scipy.special

#
# Predict next linear
#
def predict_next_linear(f, c, zero, distkm, cref):
    h = 0.001

    if ((f - h) < numpy.min(cref.x) or (f + h) > numpy.max(cref.x)):
        return -1.0, 0.0
    
    dcdf = (cref(f + h) - cref(f - h))/(2.0*h)

    fnext = ((c - f*dcdf)*zero)/(2.0*numpy.pi*distkm - dcdf*zero)
    cnext = 2.0*numpy.pi*fnext*distkm/zero
    
    cnextap = c + (fnext - f)*dcdf

    return fnext, cnext

def predict_next(f, c, zero, distkm, cref):
    h = 0.001

    if ((f - h) < numpy.min(cref.x) or (f + h) > numpy.max(cref.x)):
        return -1.0, 0.0
    
    flin, clin = predict_next_linear(f, c, zero, distkm, cref)

    dc = c - cref(f)
    if (flin < cref.x[0] or flin > cref.x[-1]):
        return -1.0, 0.0
    
    c2 = cref(flin) + dc
    
    dcdf = (cref(f + h) - cref(f - h))/(2.0*h)

    if (numpy.abs(f - flin) < 1.0e-9):
        return flin, clin
    
    A = numpy.array([[f*f, f, 1.0],
                     [flin*flin, flin, 1.0],
                     [2.0*f, 1.0, 0]])
    b = numpy.array([c, c2, dcdf])
    q = numpy.linalg.solve(A, b)

    factor = zero/(2.0*numpy.pi*distkm)
    qa = factor * q[0]
    qb = q[1]*factor - 1.0
    qc = q[2]*factor

    qd = qb*qb - 4.0*qa*qc
    if qd < 0.0:
        print(qd, qa, qb, qc, f, c, flin, clin, c2, dcdf)
        raise Exception('No solution')

    qd = numpy.sqrt(qd)
    f1 = (-qb - qd)/(2.0*qa)
    f2 = (-qb + qd)/(2.0*qa)

    if (flin < f):
        if (numpy.abs(f1 - f) < 1.0e-9):
            fnext = f1
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (numpy.abs(f2 - f) < 1.0e-9):
            fnext = f2
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f1 < f and f2 > f):
            fnext = f1
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f2 < f and f1 > f):
            fnext = f2
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f1 < f and f2 < f):
            if (f1 > f2):
                fnext = f1
                cnext = 2.0*numpy.pi*fnext*distkm/zero
            else:
                fnext = f2
                cnext = 2.0*numpy.pi*fnext*distkm/zero
        else:
            raise Exception('Unhandled %f : %f %f' % (f, f1, f2))
    else:
        if (numpy.abs(f1 - f) < 1.0e-9):
            fnext = f1
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (numpy.abs(f2 - f) < 1.0e-9):
            fnext = f2
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f1 < f and f2 > f):
            fnext = f2
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f2 < f and f1 > f):
            fnext = f1
            cnext = 2.0*numpy.pi*fnext*distkm/zero
        elif (f1 > f and f2 > f):
            if (f1 < f2):
                fnext = f1
                cnext = 2.0*numpy.pi*fnext*distkm/zero
            else:
                fnext = f2
                cnext = 2.0*numpy.pi*fnext*distkm/zero
        else:
            raise Exception('Unhandled %.15e : %.15e %f' % (f, f1, f2))

    return fnext, cnext
        
        
    

#
# Given peaks/troughs and offset, determine the next peak/trough starting
# at i in batches. If it looks like a peak/trough has been missed, add one
# to fill in. If it looks like a spurious peak/trough has been incorrectly
# identified, remove it/them.
#
def fix_forward_step(cref, distkm, j1zeros, offset, i, batches, fixed_batches):

    #
    # Compute the preceding two peaks/troughs assuming these are
    # "correct".
    #
    j = len(fixed_batches)

    fd2 = numpy.mean(fixed_batches[-2][1])
    fd1 = numpy.mean(fixed_batches[-1][1])

    cd2 = 2.0*numpy.pi*fd2 * distkm/j1zeros[offset + j - 2]
    cd1 = 2.0*numpy.pi*fd1 * distkm/j1zeros[offset + j - 1]

    fp = numpy.mean(batches[i][1])
    cp = 2.0*numpy.pi*fp * distkm/j1zeros[offset + j]

    g12 = (cd1 - cd2)/(fd1 - fd2)
    gp = (cp - cd1)/(fp - fd1)
    est_nextf, est_nextc = predict_next(fd1, cd1, j1zeros[offset + j], distkm, cref)

    halfperiod = est_nextf - fd1

    #
    # If the next frequency bin is too close (spurious peak/trough)
    # check if the next peak/trough is a better fit and skip if
    # so (or there are no more data).
    #
    if (fp < (fd1 + 0.4 * halfperiod)):

        if ((i + 2) < len(batches)):
            #
            # Check if next peak/trough is better fit
            #
            fpn = numpy.mean(batches[i + 2][1])
            if (numpy.abs(fp - est_nextf) > numpy.abs(fpn - est_nextf)):
                #
                # Skip
                #
                print('Repairing too short (next better)')
                return i + 2, fixed_batches
        else:
            print('Repairing too short', fp, fd1, fd1 + 0.5 * halfperiod)
            return i + 2, fixed_batches
    #
    # If the next frequency bin is too far, this batch may be a
    # merge of two neighboring peaks/troughs with a missing
    # intermediate trough/peak
    #
    if (fp > (fd1 + 1.6 * halfperiod)):

        #
        # Check if the batch is spread across the next trough/peak
        #
        fmin = numpy.min(batches[i][1])
        fmax = numpy.max(batches[i][1])

        if (fmin < est_nextf and fmax > est_nextf):
            #
            # Split
            #
            up, fps = batches[i]
            fleft = filter(lambda x: x < est_nextf, fps)
            fright = filter(lambda x: x > est_nextf, fps)

            fixed_batches.append((up, fleft))
            fixed_batches.append((not up, [est_nextf]))
            fixed_batches.append((up, fright))

            return i + 1, fixed_batches

        else:
            if (fmin < (fd1 + 1.6 * halfperiod)):

                up, fps = batches[i]
                fleft = filter(lambda x: x < (fd1 + 1.6 * halfperiod), fps)
                fixed_batches.append((up, fleft))
                return i + 1, fixed_batches

            else:
                print('Not split but fat', fp, fmin, fmax, est_nextf, (fd1 + 1.6 * halfperiod))
        
        
        
    #
    # Normal add
    #
    fixed_batches.append(batches[i])
    return i + 1, fixed_batches

def fix_forward(f0, j1zeros, batches, distkm, phaseref):
    
    #
    # First find nearest peak to target phase velocity
    #
    best = 1.0e9
    besti = -1
    bestf = 0.0
    for i, (up, fs) in enumerate(batches):

        f = numpy.mean(fs)
        if up:
            d = numpy.abs(f - f0)
            if (d < best):
                best = d;
                besti = i
                bestf = f

    #
    # Now determine best offset 
    #
    c0 = phaseref(bestf)

    offset = 0
    x = j1zeros[offset]
    if (-scipy.special.j1(x) < 0.0):
        offset = offset + 1
    c = 2.0*numpy.pi*bestf * distkm/j1zeros[offset]
    bestd = numpy.abs(c - c0)
    while (True):
        
        tc = 2.0*numpy.pi*bestf * distkm/j1zeros[offset + 2]
        d = numpy.abs(tc - c0)

        if (d < bestd):
            bestd = d
            c = tc
            offset = offset + 2
        else:
            break

    #
    # Check if previous peak is too close
    #
    fpp = numpy.mean(batches[besti - 1][1])
    est_prevf, est_prevc = predict_next(bestf, c, j1zeros[offset - 1], distkm, phaseref)
    ratio = (bestf - fpp)/(bestf - est_prevf)
    while (ratio < 0.4):
        offset = offset + 2
        besti = besti + 2
        bestf = numpy.mean(batches[besti][1])
        c = 2.0*numpy.pi*bestf * distkm/j1zeros[offset]
        fpp = numpy.mean(batches[besti - 1][1])
        est_prevf, est_prevc = predict_next(bestf, c, j1zeros[offset - 1], distkm, phaseref)
        ratio = (bestf - fpp)/(bestf - est_prevf)
        
    #
    # Offset is offset from i, need to subtract i
    #
    offset = offset - besti

    #
    # Now check for missing peaks and troughs
    #

    i = besti
    L = len(batches)

    #
    # We assume this is a good peak and start from here
    #
    fixed_batches = batches[:i + 1]
    i = i + 1

    while i < L:
        i, fixedbatches = fix_forward_step(phaseref, distkm, j1zeros, offset, i, batches, fixed_batches)

    return fixed_batches

def fix_backward_step(cref, distkm, j1zeros, offset, i, batches, fixed_batches, fixed_batches_start):

    #
    # Compute the preceding two peaks/troughs assuming these are
    # "correct".
    #
    j = fixed_batches_start - len(fixed_batches) - 1
    if (offset + j < 0):
        return -1, fixed_batches
    
    fd2 = numpy.mean(fixed_batches[1][1])
    fd1 = numpy.mean(fixed_batches[0][1])

    cd2 = 2.0*numpy.pi*fd2 * distkm/j1zeros[offset + j + 2]
    cd1 = 2.0*numpy.pi*fd1 * distkm/j1zeros[offset + j + 1]

    fp = numpy.mean(batches[i][1])
    cp = 2.0*numpy.pi*fp * distkm/j1zeros[offset + j]

    est_nextf, est_nextc = predict_next(fd1, cd1, j1zeros[offset + j], distkm, cref)

    if (fp > fd1):
        return i - 2, fixed_batches


    halfperiod = fd1 - est_nextf

    if (fp > fd1 - halfperiod*0.4):

        print('skipping', fp, fd1, fd1 - halfperiod*0.4, i)
        if (i > 2):
            print('', numpy.mean(batches[i - 1][1]), numpy.mean(batches[i - 2][1]))
        return i - 2, fixed_batches

    elif (fp > (fd1 - halfperiod) and i > 2):

        # Check if the next peak/trough would be better
        fp2 = numpy.mean(batches[i - 2][1])
        if (numpy.abs(fp2 - est_nextf) < numpy.abs(fp - est_nextf)):
            # Skip
            print('skipping 2', fp, fd1 - halfperiod, est_nextf, fp2)
            return i - 2, fixed_batches

    elif (fp < fd1 - halfperiod*1.6):

        est_nextf2 = cd1*j1zeros[offset + j - 1]/(2.0*numpy.pi*distkm)
        if (fp < est_nextf2):
            #
            # Insert
            #
            up, _ = batches[i]
            fixed_batches.insert(0, (up, [est_nextf]))
            fixed_batches.insert(0, (not up, [est_nextf2]))
            return i, fixed_batches
        else:
            fmin = numpy.min(batches[i][1])
            fmax = numpy.max(batches[i][1])
            if (fmax > fd1 - halfperiod*1.6):
                up, fps = batches[i]
                fright = filter(lambda x: x > (fd1 - halfperiod*1.6), fps)

                fixed_batches.insert(0, (up, fright))
                return i - 1, fixed_batches

            elif (numpy.abs((fd1 - fp)/halfperiod - 3.0) < numpy.abs((fd1 - fp)/halfperiod - 1.0)):

                up, fs = batches[i]
                fixed_batches.insert(0, (up, [est_nextf]))
                fixed_batches.insert(0, (not up, [est_nextf2]))

                return i, fixed_batches

            else:
                #
                # Ignore
                #
                up, fs = batches[i]
                delta = (fd1 - fp)/3.0
                fixed_batches.insert(0, (up, [est_nextf]))
                fixed_batches.insert(0, (not up, [est_nextf2]))
                return i, fixed_batches
                
                
            print('unfixed: missing', j, fp, est_nextf, fd1 - halfperiod*1.6, fmax, (fd1 - fp)/halfperiod, est_nextf2, halfperiod)

    fixed_batches.insert(0, batches[i])
    return i - 1, fixed_batches
    
def append_backward_step(cref, distkm, j1zeros, offset, fixed_batches, fixed_batches_start):

    #
    # Compute the preceding two peaks/troughs assuming these are
    # "correct".
    #
    j = fixed_batches_start - len(fixed_batches) - 1

    fd2 = numpy.mean(fixed_batches[1][1])
    fd1 = numpy.mean(fixed_batches[0][1])

    cd2 = 2.0*numpy.pi*fd2 * distkm/j1zeros[offset + j + 2]
    cd1 = 2.0*numpy.pi*fd1 * distkm/j1zeros[offset + j + 1]

    if (offset + j > 0):
        est_nextf, est_nextc = predict_next(fd1, cd1, j1zeros[offset + j], distkm, cref)

        if est_nextf > 0.0:
            up, _ = fixed_batches[0]
            fixed_batches.insert(0, (not up, [est_nextf]))
            return True

    return False
    
def fix_backward(f0, j1zeros, batches, distkm, phaseref):
    
    #
    # First find nearest peak to target phase velocity
    #
    best = 1.0e9
    besti = -1
    bestf = 0.0
    for i, (up, fs) in enumerate(batches):

        f = numpy.mean(fs)
        if up:
            d = numpy.abs(f - f0)
            if (d < best):
                best = d;
                besti = i
                bestf = f

    #
    # Now determine best offset 
    #
    c0 = phaseref(bestf)

    offset = 0
    x = j1zeros[offset]
    if (-scipy.special.j1(x) < 0.0):
        offset = offset + 1
    c = 2.0*numpy.pi*bestf * distkm/j1zeros[offset]
    bestd = numpy.abs(c - c0)
    while (True):
        
        tc = 2.0*numpy.pi*bestf * distkm/j1zeros[offset + 2]
        d = numpy.abs(tc - c0)

        if (d < bestd):
            bestd = d
            c = tc
            offset = offset + 2
        else:
            break

    #
    # Check if previous peak is too close
    #
    fpp = numpy.mean(batches[besti - 1][1])
    est_prevf, est_prevc = predict_next(bestf, c, j1zeros[offset - 1], distkm, phaseref)
    ratio = (bestf - fpp)/(bestf - est_prevf)
    while (ratio < 0.4):
        offset = offset + 2
        besti = besti + 2
        bestf = numpy.mean(batches[besti][1])
        c = 2.0*numpy.pi*bestf * distkm/j1zeros[offset]
        fpp = numpy.mean(batches[besti - 1][1])
        est_prevf, est_prevc = predict_next(bestf, c, j1zeros[offset - 1], distkm, phaseref)
        ratio = (bestf - fpp)/(bestf - est_prevf)

    #
    # Now check for missing peaks and troughs
    #

    i = besti
    L = len(batches)

    #
    # We assume this is a good peak and start from here
    #
    fixed_batches = batches[i:]
    fixed_start_length = len(fixed_batches)
    i = i - 1
    
    while i >= 0:
        i, fixedbatches = fix_backward_step(phaseref, distkm, j1zeros, offset, i,
                                            batches, fixed_batches, fixed_start_length)

    while append_backward_step(phaseref, distkm, j1zeros, offset, fixed_batches, fixed_start_length):
        pass
    
    return fixed_batches


if __name__ == '__main__':

    import numpy
    import scipy.interpolate
    loveref = numpy.loadtxt('../Reference/reference/reference_love_fine.txt', skiprows = 1)
    loveref = scipy.interpolate.interp1d(loveref[:,0], loveref[:,1]/1.0e3)
    
    j1zeros = scipy.special.jn_zeros(1, 500)

    offset = 0
    f = 0.10
    c = loveref(f)
    distkm = 100.0

    bestd = 1.0e9
    bestf = f
    for i, z in enumerate(j1zeros):

        fz = c/(2.0*numpy.pi * distkm) * j1zeros[i]
        d = numpy.abs(fz - f)
        if (d < bestd):
            offset = i
            bestd = d
            bestf = fz

    print('Target', f, c)
    f = bestf
    c = 2.0*numpy.pi*f * distkm/j1zeros[offset]
    print('Got', f, c, loveref(f))

    fn, cn = predict_next(f, c, j1zeros[offset + 1], distkm, loveref)

    print(f, c)
    print(fn, cn, loveref(fn), fn - f)

    fn, cn = predict_next(f, c, j1zeros[offset - 1], distkm, loveref)
            
    print(fn, cn, loveref(fn), f - fn)
    
