import sys
import os

import argparse
import numpy

import scipy.signal

import scipy.special

import matplotlib.pyplot as P

import estimaterepair

MAX_GRADIENT_DEVIATION = 5.0

WINDOW_HALF_WIDTH = 0.5

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
def loaddispersion(fname):

    f = open(fname, 'r')
    lines = f.readlines()

    slon, slat, dlon, dlat, distkm = map(float, lines[0].split())
    freq, count, acsn, csn, N = map(float, lines[1].split())

    f, r, i, ncfr, ncfi = zip(*map(lambda x: map(float, x.split()), lines[2:]))

    spec = numpy.array(r) + numpy.array(i)*1.0j
    ncf = numpy.array(ncfr) + numpy.array(ncfi)*1.0j
    return (slon, slat, dlon, dlat, distkm, int(count)), numpy.array(f), freq, acsn, csn, spec, ncf

def validate(points, phaseref):

    s, f, c, o = zip(*points)

    p = numpy.poly1d(numpy.polyfit(f, c, 2))
    pd = p.deriv()

    if (numpy.max(c) > 6.0):
        print(bcolors.FAIL + '  Invalid: maxium phase: ' + bcolors.ENDC, numpy.max(c))
        return False, 2

    return True, 0

def find_reference_trough(zero, distkm, freq, cref):

    f0 = numpy.min(cref.x)
    fl = numpy.max(cref.x)
    
    i = numpy.where(freq >= f0)[0][0]

    A = zero/(2.0 * numpy.pi * distkm)
    
    res = freq[i] - cref(freq[i]) * A
    f = freq[i]
    fi = i
    i = i + 1

    while (freq[i] <= fl):

        tres = freq[i] - cref(freq[i]) * A

        if (numpy.abs(tres) < numpy.abs(res)):
            res = tres
            f = freq[i]
            fi = i
            
        i = i + 1

    return fi, f, res

def estimate_first_trough_offset(j1zeros, points, distkm, freq, signal, cref):

    maxamplitude = numpy.max(numpy.abs(signal))
    ampthreshold = maxamplitude * 0.25

    #
    # Find first trough
    #
    for i, (s, f, c, o) in enumerate(points):
        if s == -1:

            ai = numpy.abs(f - freq).argmin()
            if signal[ai] < -ampthreshold:
                break

    #
    # Extra check: if first trough is higher frequency, reduce threshold
    # to try and get one closer to f = 0
    #
    if (f > 0.10):
        ampthreshold = maxamplitude * 0.10
        #
        # Find first trough
        #
        for i, (s, f, c, o) in enumerate(points):
            if s == -1:
                
                ai = numpy.abs(f - freq).argmin()
                if signal[ai] < -ampthreshold:
                    break
        
            

    #
    # Esimate next
    #
    est_nextzf, est_nextzc = estimaterepair.predict_next(f, c, j1zeros[o + 2], distkm, cref)

    #
    #
    #
    width = est_nextzf - f
    n = numpy.floor((f - width/2.0)/width)
    tf = f - (float(n) * width)

    print('  First detected trough at: %15.9f %d' % (f, o))
    print('  First estimated trough  : %15.9f %d' % (tf, o - n*2))

    orig = float(o)/2.0 * width + width/2.0
    print('  Origin est              : %15.9f %15.9f' % (f - orig, (f - orig)/width))
    
    #
    # Estimate using reference
    #
    refc = cref(f)

    refest_nextzf, refest_nextzc = estimaterepair.predict_next(f, refc, j1zeros[o + 2], distkm, cref)
    width = est_nextzf - f
    n = numpy.floor((f - width/2.0)/width)
    tf = f - (float(n) * width)

    print('  First estimated trough r: %15.9f %d' % (tf, o - n*2))

    best_offset = o
    fi, fref, res = find_reference_trough(j1zeros[best_offset], distkm, freq, cref)
    best_dist = numpy.abs(f - fref)

    score = best_dist
    delta_offset = 0
    
    print('  Reference offset pred: %15.9f %15.9f %15.9f' % (f, fref, best_dist))

    if (f > fref):

        # Try +ve offsets
        trial_offset = best_offset + 2
        trial_dist = 1e9
        while True:
            fi, fref, res = find_reference_trough(j1zeros[trial_offset], distkm, freq, cref)
            trial_dist = numpy.abs(f - fref)

            print('    up Trial offset pred: %d %15.9f %15.9f' % (trial_offset, fref, trial_dist))
            if (trial_dist < best_dist):
                best_offset = trial_offset
                best_dist = trial_dist
                trial_offset = trial_offset + 2
                delta_offset = delta_offset + 2

            else:
                break

        if (True or trial_dist/best_dist < 2.0):
            print('    Almost between two reference peaks', delta_offset, score, delta_offset + 2, trial_dist)
            return delta_offset, score, (delta_offset, best_dist, delta_offset + 2, trial_dist)

    else:

        # Try -ve offsets
        trial_offset = best_offset - 2
        trial_dist = 1e9
        while trial_offset >= 0:
            fi, fref, res = find_reference_trough(j1zeros[trial_offset], distkm, freq, cref)
            trial_dist = numpy.abs(f - fref)

            print('    dn Trial offset pred: %d %15.9f %15.9f' % (trial_offset, fref, trial_dist))
            if (trial_dist < best_dist):
                best_offset = trial_offset
                best_dist = trial_dist
                trial_offset = trial_offset - 2
                delta_offset = delta_offset - 2

            else:
                break

        if (True or trial_dist/best_dist < 2.0):
            print('    Almost between two reference peaks', delta_offset, score, delta_offset - 2, trial_dist)
            return delta_offset, score, (delta_offset, best_dist, delta_offset - 2, trial_dist)

        
    #find_reference_trough(zero, distkm, freq, cref):
    print('  Recommended first trough offset: %d (%d) %15.9f' % (best_offset, delta_offset, score))
    return delta_offset, score, None

def lstscore(points, phaseref):

    s, f, c, o = zip(*points)
    p = numpy.poly1d(numpy.polyfit(f, c, 6))

    xmin, xmax = 0.05, 0.075

    if min(f) > xmin or max(f) < xmax:
        print('exiting')
        return 1.0e9
    
    x = numpy.linspace(xmin, xmax, 32)
    score = numpy.sqrt(numpy.sum((p(x) - phaseref(x))**2))

    print(score)

    pd = p.deriv()
    h = 0.001
    dscore = 0.0
    for ix in x:
        d = pd(ix) - (phaseref(ix + h) - phaseref(ix - h))/(2.0 * h)
        dscore = dscore + d*d

    print(dscore)
        
    return score #+ 2.0e-1*numpy.sqrt(dscore)

def findpeak(f, signal, if0, if1):
    if if0 < 0:
        if0 = 0
    if if1 >= f.size:
        if1 = f.size - 1
    if if0 >= if1:
        return -1

    mi = numpy.argmax(signal[if0:if1 + 1])
    mi = if0 + mi
    
    return mi

def findzerocross(sign, f, signal, if0, if1):
    if if0 < 0:
        if0 = 0
    if if1 >= f.size:
        if1 = f.size - 1
    if if0 >= if1:
        return -1
    
    t = signal[if0:if1]
    t2 = signal[if0 + 1: if1 + 1]
    indices = numpy.where((t * t2) < 0.0)[0]
    
    if indices.size == 0:
        # None found
        return -1

    elif indices.size == 1:
        # One found
        zci = indices[0] + if0
        if (signal[zci] < 0.0 and signal[zci+1] > 0.0):
            if sign < 0:
                return zci
            else:
                return -1
        elif (signal[zci] > 0.0 and signal[zci+1] < 0.0):
            if sign > 0:
                return zci
            else:
                return -1

        else:
            raise Exception('Invalid zero cross')

    else:
        # Assume if0, if1 represent extremal bounds and select nearest to centre
        besti = -1
        bestscore = 8192
        target = (if1 - if0)//2
        for i, ti in enumerate(indices):
            score = numpy.abs(ti - target)
            zci = ti + if0
            if sign < 0:
                if (signal[zci] < 0.0 and signal[zci+1] > 0.0):
                    if score < bestscore:
                        bestscore = score
                        besti = i
            elif sign > 0:
                if (signal[zci] > 0.0 and signal[zci+1] < 0.0):
                    if score < bestscore:
                        bestscore = score
                        besti = i

        if besti < 0:
            return -1
        else:
            zci = indices[besti] + if0
            return zci

def findtrough(f, signal, if0, if1):
    if if0 < 0:
        if0 = 0
    if if1 >= f.size:
        if1 = f.size - 1
    if if0 >= if1:
        return -1
        
    mi = numpy.argmin(signal[if0:if1 + 1])
    mi = if0 + mi
    
    return mi

def mkwindow(freq, fmin, fmax):

    if0 = numpy.where(freq > fmin)[0][0] - 1
    if1 = numpy.where(freq > fmax)[0][0] + 1

    return if0, if1

################################################################################
#
# Forward
#
################################################################################

#
# Recursively try to find next peak forward, followed by next trough if unable to
# find peak, ignoring zero crossings.
#
def find_forward_peak(j0zeros, j1zeros, f, c, offset,
                      freq, signal, distkm, maxamplitude, cref,
                      fmin, fmax, threshold):

    if offset % 2 != 1:
        raise Exception('Look for peak with even offset')
    
    #
    # Estimate next peak location
    #
    est_nextpf, est_nextpc = estimaterepair.predict_next(f, c, j1zeros[offset], distkm, cref)
    if est_nextpf < 0.0:
        return False, 0, 0.0, 0.0, -1
    
    #
    # Estimate next zero cross to get window width (ie, peak should be in range est_f +- zero cross)
    #
    est_nextzf, est_nextzc = estimaterepair.predict_next(f, c, j0zeros[offset], distkm, cref)

    if est_nextzf > est_nextpf:
        print(f, c, offset)
        print(est_nextpf, est_nextpc)
        print(est_nextzf, est_nextzc)
        raise Exception('Unexpected ordering of peak/zero')
    
    wfwidth = est_nextpf - est_nextzf 
    wfmin = est_nextpf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextpf + wfwidth*WINDOW_HALF_WIDTH
    if0, if1 = mkwindow(freq, wfmin, wfmax)
    pci = findpeak(freq, signal, if0, if1)

    if pci > 0 and (signal[pci] > maxamplitude*threshold):

        if (pci == if0):
            # At edge, check further
            tpci = pci
            while (tpci > 0 and signal[tpci - 1] > signal[tpci]):
                tpci = tpci - 1

            binthreshold = (if1 - if0)//2
            if (True or pci - tpci < binthreshold):
                #print 'Warning: peak slightly out of range lower', pci, tpci
                pci = tpci

        elif (pci == if1):
            tpci = pci
            while (tpci < freq.size - 1 and signal[tpci + 1] > signal[tpci]):
                tpci = tpci + 1

            binthreshold = (if1 - if0)//2
            if (True or tpci - pci < binthreshold):
                #print 'Warning: peak slightly out of range higher', pci, tpci
                pci = tpci

        #
        # Found
        #
        next_f = freq[pci]
        next_c = 2.0*numpy.pi*next_f * distkm/j1zeros[offset]

        if (next_f > fmax):
            # Out of range
            print(bcolors.WARNING + '  Ignoring peak: out of range' + bcolors.ENDC)
            return False, 1, next_f, next_c, offset

        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c > MAX_GRADIENT_DEVIATION):
            # Gradient deviation
            pci = -1

        else:
            return True, 1, next_f, next_c, offset

    #
    # Try for next trough
    #
    return find_forward_trough(j0zeros, j1zeros, f, c, offset + 1,
                               freq, signal, distkm, maxamplitude, cref,
                               fmin, fmax, threshold)

#
# Recursively try to find next peak forward, followed by next trough if unable to
# find peak, ignoring zero crossings.
#
def find_forward_trough(j0zeros, j1zeros, f, c, offset,
                        freq, signal, distkm, maxamplitude, cref,
                        fmin, fmax, threshold):

    if offset % 2 != 0:
        raise Exception('Look for trough with odd offset')
    
    #
    # Estimate next peak location
    #
    est_nexttf, est_nexttc = estimaterepair.predict_next(f, c, j1zeros[offset], distkm, cref)
    if est_nexttf < 0.0:
        return False, 0, 0.0, 0.0, -1
    
    #
    # Estimate next zero cross to get window width (ie, peak should be in range est_f +- zero cross)
    #
    est_nextzf, est_nextzc = estimaterepair.predict_next(f, c, j0zeros[offset], distkm, cref)

    if est_nextzf > est_nexttf:
        raise Exception('Unexpected ordering of trough/zero')
    
    wfwidth = est_nexttf - est_nextzf 
    wfmin = est_nexttf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nexttf + wfwidth*WINDOW_HALF_WIDTH
    if0, if1 = mkwindow(freq, wfmin, wfmax)
    pci = findtrough(freq, signal, if0, if1)
                
    if pci > 0 and (signal[pci] < -maxamplitude*threshold):

        if (pci == if0):
            # At edge, check further
            tpci = pci
            while (tpci > 0 and signal[tpci - 1] < signal[tpci]):
                tpci = tpci - 1

            binthreshold = (if1 - if0)//2
            if (True or pci - tpci < binthreshold):
                print('Warning: trough slightly out of range lower', pci, tpci)
                pci = tpci

        elif (pci == if1):
            tpci = pci
            while (tpci < freq.size - 1 and signal[tpci + 1] < signal[tpci]):
                tpci = tpci + 1

            binthreshold = (if1 - if0)//2
            if (True or tpci - pci < binthreshold):
                print('Warning: trough slightly out of range higher', pci, tpci)
                pci = tpci
                

                
        #
        # Found
        #
        next_f = freq[pci]
        next_c = 2.0*numpy.pi*next_f * distkm/j1zeros[offset]

        if (next_f > fmax):
            # Out of range
            return False, -1, next_f, next_c, offset

        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c > MAX_GRADIENT_DEVIATION):
            # Gradient deviation
            pci = -1

        else:
            return True, -1, next_f, next_c, offset

    #
    # Try for next peak
    #
    return find_forward_peak(j0zeros, j1zeros, f, c, offset + 1,
                             freq, signal, distkm, maxamplitude, cref,
                             fmin, fmax, threshold)

def add_next_forward_from_peak(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref,
                               fmin, fmax,
                               picks,
                               threshold):
    
    #
    # Starting from peak the algorithm is:
    #
    # 1. Attempt to find downward (+ve - -ve) zero cross, if ok use, otherwise
    # 2. Attempt to find next trough, if ok use, otherwise
    # 3. Attempt to find next peak
    #
    # Reasons for rejection:
    #
    # 1. Sharp decrease in phase velocity -> spurious peak/trough/cross
    # 2. Sharp increase in phase velocity -> missing peak/trough/cross
    # 3. Low magnitude peak/trough
    #
    
    sign, f, c, offset = picks[-1]

    next_offset = offset + 1 # Offset for next zero cross and trough
    est_nextf, est_nextc = estimaterepair.predict_next(f, c, j0zeros[next_offset], distkm, cref)

    if est_nextf < f:
        #
        # Can happen near end
        return picks, True

    if est_nextf > fmax:
        return picks, True

    # Window for crossing
    wfwidth = est_nextf - f
    wfmin = est_nextf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextf + wfwidth*WINDOW_HALF_WIDTH
        
    if0, if1 = mkwindow(freq, wfmin, wfmax)

    zci = findzerocross(sign, freq, signal, if0, if1)
    if (zci > 0):
        t = -signal[zci]/(signal[zci + 1] - signal[zci])
        next_f = freq[zci] + (freq[zci + 1] - freq[zci])*t
        next_c = 2.0*numpy.pi*next_f * distkm/j0zeros[next_offset]

        deltac = next_c - c
        if f > 0.1 and numpy.abs(deltac) > 0.1:
            # Strong positive/negative change in phase, invalidate and fall through
            zci = -1
        else:
            picks.append((0, next_f, next_c, next_offset))
            return picks, False

    #
    # Fall through: try to find next trough and recursively next peak etc
    #
    found, up, ff, fc, foffset = find_forward_trough(j0zeros, j1zeros, f, c, next_offset,
                                                     freq, signal, distkm, maxamplitude, cref,
                                                     fmin, fmax, threshold)
    
    if found:
        picks.append((up, ff, fc, foffset))
        return picks, False

    return picks, True
        
def add_next_forward_from_trough(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref,
                                 fmin, fmax,
                                 picks,
                                 threshold):
    
    #
    # Starting from trough the algorithm is:
    #
    # 1. Attempt to find upward (-ve - +ve) zero cross, if ok use, otherwise
    # 2. Attempt to find next peak, if ok use, otherwise
    # 3. Attempt to find next trough, if ok use, otherwise
    # 4. Goto 2 util fmax
    #
    # Reasons for rejection:
    #
    # 1. Sharp decrease in phase velocity -> spurious peak/trough/cross
    # 2. Sharp increase in phase velocity -> missing peak/trough/cross
    # 3. Low magnitude peak/trough
    #
    
    sign, f, c, offset = picks[-1]

    next_offset = offset + 1 # Offset for next zero cross and trough
    est_nextf, est_nextc = estimaterepair.predict_next(f, c, j0zeros[next_offset], distkm, cref)

    if est_nextf < f:
        #
        # Can happen near end
        return picks, True

    if est_nextf > fmax:
        return picks, True

    # Window for crossing
    wfwidth = est_nextf - f
    wfmin = est_nextf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextf + wfwidth*WINDOW_HALF_WIDTH
        
    if0, if1 = mkwindow(freq, wfmin, wfmax)

    zci = findzerocross(sign, freq, signal, if0, if1)
    if (zci > 0):
        t = -signal[zci]/(signal[zci + 1] - signal[zci])
        next_f = freq[zci] + (freq[zci + 1] - freq[zci])*t
        next_c = 2.0*numpy.pi*next_f * distkm/j0zeros[next_offset]

        deltac = next_c - c
        if f > 0.1 and numpy.abs(deltac) > 0.1:
            # Strong positive/negative change in phase, invalidate and fall through
            zci = -1
        else:
            picks.append((0, next_f, next_c, next_offset))
            return picks, False

    #
    # Fall through: try to find next trough and recursively next peak etc
    #
    found, up, ff, fc, foffset = find_forward_peak(j0zeros, j1zeros, f, c, next_offset,
                                                   freq, signal, distkm, maxamplitude, cref,
                                                   fmin, fmax, threshold)
    
    if found:
        picks.append((up, ff, fc, foffset))
        return picks, False

    # Not found
    return picks, True
    
def add_next_forward(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref, fmin, fmax, picks, threshold):

    sign, f, c, offset = picks[-1]

    if sign == 1: # Expect nve zero cross

        return add_next_forward_from_peak(j0zeros, j1zeros,
                                          freq, signal, distkm, maxamplitude, cref,
                                          fmin, fmax,
                                          picks,
                                          threshold)
        

    elif sign == 0: # Expect trough if offset even, peak if odd

        if offset % 2 == 1:

            # Zero cross -> peak, offset is unchanged
            found, up, ff, fc, foffset = find_forward_peak(j0zeros, j1zeros, f, c, offset,
                                                           freq, signal, distkm, maxamplitude, cref,
                                                           fmin, fmax, threshold)

            if (found):
                picks.append((up, ff, fc, foffset))
                return picks, False
                
        else:
            # Zero cross -> trough, offset is unchanged
            found, up, ff, fc, foffset = find_forward_trough(j0zeros, j1zeros, f, c, offset,
                                                             freq, signal, distkm, maxamplitude, cref,
                                                             fmin, fmax, threshold)

            if (found):
                picks.append((up, ff, fc, foffset))
                return picks, False
    
    else: # Trough -> Expect pve zero cross

        return add_next_forward_from_trough(j0zeros, j1zeros,
                                            freq, signal, distkm, maxamplitude, cref,
                                            fmin, fmax,
                                            picks,
                                            threshold)
    

    return picks, True

################################################################################
#
# Backward
#
################################################################################

#
# Recursively try to find previous peak backward, followed by previous trough if unable to
# find peak, ignoring zero crossings.
#
def find_backward_peak(j0zeros, j1zeros, f, c, offset,
                       freq, signal, distkm, maxamplitude, cref,
                       fmin, fmax, threshold):
    
    if offset % 2 != 1:
        raise Exception('Look for peak with even offset')
    
    #
    # Estimate previous peak location
    #
    est_nextpf, est_nextpc = estimaterepair.predict_next(f, c, j1zeros[offset], distkm, cref)
    if est_nextpf < 0.0:
        return False, 0, 0.0, 0.0, -1
    
    #
    # Estimate next zero cross to get window width (ie, peak should be in range est_f +- zero cross)
    #
    est_nextzf, est_nextzc = estimaterepair.predict_next(f, c, j0zeros[offset + 1], distkm, cref)

    if est_nextpf > est_nextzf:
        print(f, c, offset)
        print(est_nextpf, est_nextpc)
        print(est_nextzf, est_nextzc)
        raise Exception('Unexpected ordering of peak/zero')
    
    wfwidth = est_nextzf - est_nextpf 
    wfmin = est_nextpf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextpf + wfwidth*WINDOW_HALF_WIDTH
    if0, if1 = mkwindow(freq, wfmin, wfmax)
    pci = findpeak(freq, signal, if0, if1)

    if pci > 0 and (signal[pci] > maxamplitude*threshold):

        if (pci == if0):
            # At edge, check further
            tpci = pci
            while (tpci > 0 and signal[tpci - 1] > signal[tpci]):
                tpci = tpci - 1

            binthreshold = (if1 - if0)//2
            if (pci - tpci < binthreshold):
                #print '  find_backward_peak: Warning: peak slightly out of range lower', pci, tpci
                pci = tpci

        elif (pci == if1):
            tpci = pci
            while (tpci < freq.size - 1 and signal[tpci + 1] > signal[tpci]):
                tpci = tpci + 1

            binthreshold = (if1 - if0)//2
            if (tpci - pci < binthreshold):
                #print '  find_backward_peak: Warning: peak slightly out of range higher', pci, tpci
                pci = tpci

        #
        # Found
        #
        next_f = freq[pci]
        next_c = 2.0*numpy.pi*next_f * distkm/j1zeros[offset]

        if (next_f < fmin):
            # Out of range
            #print bcolors.WARNING + '  find_backward_peak: Ignoring peak: out of range' + bcolors.ENDC
            return False, 1, next_f, next_c, offset

        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c > MAX_GRADIENT_DEVIATION):
            # Gradient deviation
            #print bcolors.WARNING + '  find_backward_peak: Ignoring peak: gradient deviation: ' + bcolors.ENDC, delta_c, est_delta_c, rel_delta_c
            pci = -1

        else:
            #print '  find_backward_peak: Found peak: predicted %15.9f found %15.9f err %12.4e' % (est_nextpf, next_f,
                                                                                                  #est_nextpf - next_f)
            return True, 1, next_f, next_c, offset

    #
    # Try for next trough
    #
    return find_backward_trough(j0zeros, j1zeros, f, c, offset - 1,
                               freq, signal, distkm, maxamplitude, cref,
                               fmin, fmax, threshold)

#
# Recursively try to find previous peak backward, followed by previous trough if unable to
# find peak, ignoring zero crossings.
#
def find_backward_trough(j0zeros, j1zeros, f, c, offset,
                         freq, signal, distkm, maxamplitude, cref,
                         fmin, fmax, threshold):

    if offset % 2 != 0:
        raise Exception('Look for trough with odd offset')
    
    #
    # Estimate next peak location
    #
    est_nexttf, est_nexttc = estimaterepair.predict_next(f, c, j1zeros[offset], distkm, cref)
    if est_nexttf < 0.0:
        return False, 0, 0.0, 0.0, -1
    
    #
    # Estimate next zero cross to get window width (ie, peak should be in range est_f +- zero cross)
    #
    est_nextzf, est_nextzc = estimaterepair.predict_next(f, c, j0zeros[offset + 1], distkm, cref)

    if est_nexttf > est_nextzf:
        raise Exception('Unexpected ordering of trough/zero')
    
    wfwidth = est_nextzf - est_nexttf 
    wfmin = est_nexttf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nexttf + wfwidth*WINDOW_HALF_WIDTH
    if0, if1 = mkwindow(freq, wfmin, wfmax)
    pci = findtrough(freq, signal, if0, if1)
                
    if pci > 0 and (signal[pci] < -maxamplitude*threshold):

        if (pci == if0):
            # At edge, check further
            tpci = pci
            while (tpci > 0 and signal[tpci - 1] < signal[tpci]):
                tpci = tpci - 1

            binthreshold = (if1 - if0)//2
            if (True or pci - tpci < binthreshold):
                #print '  find_backward_trough: Warning: trough slightly out of range lower', pci, tpci
                pci = tpci

        elif (pci == if1):
            tpci = pci
            while (tpci < freq.size - 1 and signal[tpci + 1] < signal[tpci]):
                tpci = tpci + 1

            binthreshold = (if1 - if0)//2
            if (True or tpci - pci < binthreshold):
                #print '  find_backward_trough: Warning: trough slightly out of range higher', pci, tpci
                pci = tpci
                
        #
        # Found
        #
        next_f = freq[pci]
        next_c = 2.0*numpy.pi*next_f * distkm/j1zeros[offset]

        if (next_f < fmin):
            # Out of range
            #print '  find_backward_trough: Ignoring trough: out of range'
            return False, -1, next_f, next_c, offset

        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c > MAX_GRADIENT_DEVIATION):
            # Gradient deviation
            #print '  find_backward_trough: Ignoring trough: gradient deviation: ', delta_c, est_delta_c, rel_delta_c
            pci = -1

        else:
            return True, -1, next_f, next_c, offset

    #
    # Try for next peak
    #
    if offset >= 2:
        
        return find_backward_peak(j0zeros, j1zeros, f, c, offset - 1,
                                 freq, signal, distkm, maxamplitude, cref,
                                 fmin, fmax, threshold)
    else:

        return False, -1, 0.0, 0.0, 0

def add_next_backward_from_peak(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref,
                                fmin, fmax,
                                picks,
                                threshold):
    
    #
    # Starting from peak the algorithm is:
    #
    # 1. Attempt to find upward (-ve - +ve) zero cross, if ok use, otherwise
    # 2. Attempt to find previous trough, if ok use, otherwise
    # 3. Attempt to find previous peak
    #
    # Reasons for rejection:
    #
    # 1. Sharp decrease in phase velocity -> spurious peak/trough/cross
    # 2. Sharp increase in phase velocity -> missing peak/trough/cross
    # 3. Low magnitude peak/trough
    #
    
    sign, f, c, offset = picks[0]

    next_offset = offset # Offset for previous zero
    est_nextf, est_nextc = estimaterepair.predict_next(f, c, j0zeros[next_offset], distkm, cref)

    if est_nextf > f:
        #
        # Can happen near end
        return picks, True

    if est_nextf < fmin:
        return picks, True

    # Window for crossing
    wfwidth = f - est_nextf
    wfmin = est_nextf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextf + wfwidth*WINDOW_HALF_WIDTH
        
    if0, if1 = mkwindow(freq, wfmin, wfmax)

    zci = findzerocross(-sign, freq, signal, if0, if1)
    if (zci > 0):
        t = -signal[zci]/(signal[zci + 1] - signal[zci])
        next_f = freq[zci] + (freq[zci + 1] - freq[zci])*t
        next_c = 2.0*numpy.pi*next_f * distkm/j0zeros[next_offset]

        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c < MAX_GRADIENT_DEVIATION):
            picks.insert(0, (0, next_f, next_c, next_offset))
            return picks, False

    #
    # Fall through: try to find next trough and recursively next peak etc
    #
    found, up, ff, fc, foffset = find_backward_trough(j0zeros, j1zeros, f, c, next_offset - 1,
                                                      freq, signal, distkm, maxamplitude, cref,
                                                      fmin, fmax, threshold)
    if found:
        picks.insert(0, (up, ff, fc, foffset))
        return picks, False

    return picks, True
        
def add_next_backward_from_trough(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref,
                                  fmin, fmax,
                                  picks,
                                  threshold):
    
    #
    # Starting from trough the algorithm is:
    #
    # 1. Attempt to find downward (+ve - -ve) zero cross, if ok use, otherwise
    # 2. Attempt to find previous peak, if ok use, otherwise
    # 3. Attempt to find previous trough, if ok use, otherwise
    # 4. Goto 2 util fmax
    #
    # Reasons for rejection:
    #
    # 1. Sharp decrease in phase velocity -> spurious peak/trough/cross
    # 2. Sharp increase in phase velocity -> missing peak/trough/cross
    # 3. Low magnitude peak/trough
    #
    
    sign, f, c, offset = picks[0]

    next_offset = offset # Offset for next zero cross 
    est_nextf, est_nextc = estimaterepair.predict_next(f, c, j0zeros[next_offset], distkm, cref)

    if est_nextf > f:
        #
        # Can happen near end
        return picks, True

    if est_nextf < fmin:
        return picks, True

    # Window for crossing
    wfwidth = f - est_nextf 
    wfmin = est_nextf - wfwidth*WINDOW_HALF_WIDTH
    wfmax = est_nextf + wfwidth*WINDOW_HALF_WIDTH
        
    if0, if1 = mkwindow(freq, wfmin, wfmax)

    zci = findzerocross(sign, freq, signal, if0, if1)
    if (zci > 0):
        t = -signal[zci]/(signal[zci + 1] - signal[zci])
        next_f = freq[zci] + (freq[zci + 1] - freq[zci])*t
        next_c = 2.0*numpy.pi*next_f * distkm/j0zeros[next_offset]

        deltac = next_c - c
        delta_c = next_c - c
        est_delta_c = cref(next_f) - cref(f)
        rel_delta_c = numpy.abs(delta_c - est_delta_c)/numpy.abs(est_delta_c)
        if (rel_delta_c < MAX_GRADIENT_DEVIATION):
            picks.insert(0, (0, next_f, next_c, next_offset))
            return picks, False

    #
    # Fall through: try to find next trough and recursively next peak etc
    #
    if (next_offset >= 2) :
        found, up, ff, fc, foffset = find_backward_peak(j0zeros, j1zeros, f, c, next_offset - 1,
                                                        freq, signal, distkm, maxamplitude, cref,
                                                        fmin, fmax, threshold)
    
        if found:
            picks.insert(0, (up, ff, fc, foffset))
            return picks, False

    # Not found
    return picks, True



def add_next_backward(j0zeros, j1zeros, freq, signal, distkm, maxamplitude, cref, fmin, fmax, picks, threshold):

    sign, f, c, offset = picks[0]

    if sign == 1: # Expect pve zero cross

        return add_next_backward_from_peak(j0zeros, j1zeros,
                                           freq, signal, distkm, maxamplitude, cref,
                                           fmin, fmax,
                                           picks,
                                           threshold)

    elif sign == 0: # Expect trough if offset even, peak if odd

        if offset == 0: # No more zeros
            return picks, True
        
        if offset % 2 == 1:

            found, up, ff, fc, foffset = find_backward_trough(j0zeros, j1zeros,
                                                              f, c, offset - 1,
                                                              freq, signal, distkm, maxamplitude, cref,
                                                              fmin, fmax, threshold)

            if (found):
                picks.insert(0, (up, ff, fc, foffset))
                return picks, False
                
        else:
            found, up, ff, fc, foffset = find_backward_peak(j0zeros, j1zeros,
                                                            f, c, offset - 1,
                                                            freq, signal, distkm, maxamplitude, cref,
                                                            fmin, fmax, threshold)

            if (found):
                picks.insert(0, (up, ff, fc, foffset))
                return picks, False
    
    else: # Trough -> Expect nve zero cross

        return add_next_backward_from_trough(j0zeros, j1zeros,
                                             freq, signal, distkm, maxamplitude, cref,
                                             fmin, fmax,
                                             picks,
                                             threshold)


    return picks, True

def pick(j0zeros, j1zeros, freq, signal, distkm, phaseref, fmin, fmax, suggestoffset = 0, threshold = 0.075):

    #
    # Find the first peak/trough from the maximum (ie most likely to be a true peak
    #
    indices = numpy.where((freq >= 0.075) & (freq <= 0.2))[0]
    
    lp = numpy.argmax(signal[indices]) + indices[0]
    lt = numpy.argmin(signal[indices]) + indices[0]

    if (-signal[lt] > signal[lp]):

        maxamplitude = -signal[lt]
        # Trough first (troughs are even zeros of j1)
        
        f = freq[lt]

        offset = 0
        bestoffset = -1
        bestdist = 1.0e9
        cref = phaseref(f)
        #print cref
        
        while offset < len(j1zeros):
            c = 2.0*numpy.pi*f * distkm/j1zeros[offset]
            dist = numpy.abs(c - cref)
            if (dist < bestdist):
                bestdist = dist
                bestoffset = offset
                
            offset = offset + 2

        offset = bestoffset + suggestoffset
        c = 2.0*numpy.pi*f * distkm/j1zeros[offset]
        picks = [(-1.0, f, c, offset)]
        
    else:
        # Peak first (peaks are odd zeros of j1)

        maxamplitude = signal[lp]
        f = freq[lp]

        offset = 1
        bestoffset = -1
        bestdist = 1.0e9
        cref = phaseref(f)
        while offset < len(j1zeros):
            c = 2.0*numpy.pi*f * distkm/j1zeros[offset]
            dist = numpy.abs(c - cref)
            if (dist < bestdist):
                bestdist = dist
                bestoffset = offset
                
            offset = offset + 2

        offset = bestoffset + suggestoffset
        c = 2.0*numpy.pi*f * distkm/j1zeros[offset]
        picks = [(1.0, f, c, offset)]

    while True:
        picks, finished = add_next_backward(j0zeros, j1zeros,
                                            freq, signal, distkm,
                                            maxamplitude, phaseref,
                                            fmin, fmax, picks,
                                            threshold)

        if finished:
            break

    while True:

        picks, finished = add_next_forward(j0zeros, j1zeros,
                                           freq, signal, distkm,
                                           maxamplitude, phaseref,
                                           fmin, fmax, picks,
                                           threshold)

        if finished:
            break
        
    return picks

def estimate_error(j0zeros, j1zeros, s, f, c, o, distkm):

    if s == 0:

        c0 = 2.0*numpy.pi*f * distkm/j0zeros[o]
        if o >= 2:
            
            cplus = 2.0*numpy.pi*f * distkm/j0zeros[o - 2]
            err = (cplus - c0)/2.0
            
        else:

            cplus = 2.0*numpy.pi*f * distkm/j0zeros[o + 2]
            err = (cplus - c0)/2.0

    else:

        c0 = 2.0*numpy.pi*f * distkm/j1zeros[o]
        if o >= 2:
            
            cplus = 2.0*numpy.pi*f * distkm/j1zeros[o - 2]
            err = (cplus - c0)/2.0
            
        else:

            cplus = 2.0*numpy.pi*f * distkm/j1zeros[o + 2]
            err = (cplus - c0)/2.0

    return err
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--path', type = str, default = '/home/rhys/PhD/PhDS/tools/ANTEstimate/Processing', help = 'Data base path')
    parser.add_argument('-s', '--station-pair', type = str, required = True, help = 'Input station pair')

    parser.add_argument('-t', '--threshold', type = float, default = 0.05, help = 'Threshold envelope')

    
    parser.add_argument('-r', '--love-reference', type = str, default = os.path.join(sys.path[0], '../../tutorial/reference/reference_love_fine.txt'), help = 'Reference Love phase')
    parser.add_argument('-R', '--rayleigh-reference', type = str, default = os.path.join(sys.path[0], '../../tutorial/reference/reference_rayleigh_fine.txt'), help = 'Reference Rayleigh phase')

    parser.add_argument('-v', '--vmin', type = float, default = 0.5, help = 'Minimum velocity')
    parser.add_argument('-V', '--vmax', type = float, default = 7.0, help = 'Maximum velocity')

    parser.add_argument('-f', '--freq-min', type = float, default = 1.0/40.0, help = 'Min frequency')
    parser.add_argument('-F', '--freq-max', type = float, default = 0.35, help = 'Max frequency')

    parser.add_argument('-o', '--output', type = str, default = None, help = 'Output phase velocities')

    parser.add_argument('-b', '--fix-forward', action = 'store_true', default = False, help = 'Fix batches forward')
    parser.add_argument('-B', '--fix-backward', action = 'store_true', default = False, help = 'Fix batches backward')

    parser.add_argument('-O', '--offset', type = int, default = 0, help = 'Plot offset curves')

    parser.add_argument('--filter', type = float, default = 3, help = 'Filter width')

    parser.add_argument('--noshow', action = 'store_true', default = False, help = 'No plotting')
                        
    args = parser.parse_args()

    plotting = not args.noshow
    
    #
    # Load observed spectra
    #
    lovedata = os.path.join(args.path, 'LoveResponse/dispersion_%s.txt' % args.station_pair)
    (_, _, _, _, distkm, _), freq, sample_rate, loveacsn, lovecsn, lovespec, lovencf = loaddispersion(lovedata)

    rayleighdata = os.path.join(args.path, 'RayleighResponse/dispersion_%s.txt' % args.station_pair)
    (_, _, _, _, distkm, _), freq, sample_rate, rayleighacsn, rayleighcsn, rayleighspec, rayleighncf = loaddispersion(rayleighdata)
    #
    # Load reference models
    #
    loveref = numpy.loadtxt(args.love_reference, skiprows = 1)
    indices = numpy.where(loveref[:,1] > 0.0)[0]
    lovephaseref = scipy.interpolate.interp1d(loveref[indices,0], loveref[indices,1]/1.0e3)
    rayleighref = numpy.loadtxt(args.rayleigh_reference, skiprows = 1)
    indices = numpy.where(rayleighref[:,1] > 0.0)[0]
    rayleighphaseref = scipy.interpolate.interp1d(rayleighref[indices,0], rayleighref[indices,1]/1.0e3)

    #
    # Compute signals
    #
    lovesignal = numpy.real(lovencf)
    rayleighsignal = numpy.real(rayleighncf)

    if args.filter > 0.0:
        lovesignal = scipy.ndimage.gaussian_filter1d(lovesignal, args.filter)
        rayleighsignal = scipy.ndimage.gaussian_filter1d(rayleighsignal, args.filter)

    if plotting:
        fig, ax = P.subplots(2, 1)
        fig.set_tight_layout(True)
        ax[0].set_title('Love')
        ax[0].plot(freq, lovesignal)
        ax[1].set_title('Rayleigh')
        ax[1].plot(freq, rayleighsignal)
        
        ax[0].set_xlim(0, 0.5)
        ax[1].set_xlim(0, 0.5)

        fig, bx = P.subplots()
        fig.set_tight_layout(True)
        bx.set_xlim(0, 0.5)
        bx.set_ylim(0, 6)
        

    j0zeros = scipy.special.jn_zeros(0, 1024)
    j1zeros = scipy.special.jn_zeros(1, 1024)

    print('Picking Love')
    lovedoffset = 0
    alreadytried = {}
    acceptable = False
    
    while not acceptable:
        acceptable = True
        print('Love Begin pick: %d' % lovedoffset)
        lovepoints = pick(j0zeros, j1zeros, freq, lovesignal, distkm, lovephaseref,
                          args.freq_min, args.freq_max, lovedoffset)


        offset, score, lovebounds = estimate_first_trough_offset(j1zeros,
                                                                 lovepoints,
                                                                 distkm,
                                                                 freq,
                                                                 lovesignal,
                                                                 lovephaseref)
        alreadytried[lovedoffset] = (lovepoints, score)

        if offset != 0:
            acceptable = False
            print('Retrying', offset)
            lovedoffset = lovedoffset + offset
            if alreadytried.has_key(lovedoffset):
                print('Looped back on self, using best score')
                minv = 1e30
                minpts = None
                for k, (pts, v) in alreadytried.items():
                    if (v < minv):
                        minpts = pts
                        minv = v
                lovepoints = minpts
                break
        
    print('Picking Rayleigh')
    rayleighdoffset = 0
    acceptable = False
    alreadytried = {}
    
    while not acceptable:
        acceptable = True
        rayleighpoints = pick(j0zeros, j1zeros, freq, rayleighsignal, distkm, rayleighphaseref,
                              args.freq_min, args.freq_max, rayleighdoffset)

        offset, score, rayleighbounds = estimate_first_trough_offset(j1zeros,
                                                                     rayleighpoints,
                                                                     distkm,
                                                                     freq,
                                                                     rayleighsignal,
                                                                     rayleighphaseref)
        alreadytried[rayleighdoffset] = (rayleighpoints, score)

        acceptable = (offset == 0)

        if offset != 0:
            acceptable = False
            print('Retrying', offset)
            rayleighdoffset = rayleighdoffset + offset
            if alreadytried.has_key(rayleighdoffset):
                print('Looped back on self, using best score')
                minv = 1e30
                minpts = None
                for k, (pts, v) in alreadytried.items():
                    if (v < minv):
                        minpts = pts
                        minv = v
                rayleighpoints = minpts
                break


    if lovebounds is None:
        if not rayleighbounds is None:
            #
            # Special case where trough is almost halfway between reference troughs
            #
            offset1, score1, offset2, score2 = rayleighbounds
            print('Rayleigh between troughs:', rayleighdoffset, offset1, score1, offset2, score2)
            
            if offset1 == 0:
                points1 = list(rayleighpoints)
            else:
                points1 = pick(j0zeros, j1zeros, freq, rayleighsignal, distkm, rayleighphaseref,
                               args.freq_min, args.freq_max, rayleighdoffset + offset1)
            
            if offset2 == 0:
                points2 = list(rayleighpoints)
            else:
                points2 = pick(j0zeros, j1zeros, freq, rayleighsignal, distkm, rayleighphaseref,
                               args.freq_min, args.freq_max, rayleighdoffset + offset2)

            _, lf, lc, _ = zip(*lovepoints)
            lcurve = scipy.interpolate.interp1d(lf, lc)
            l10 = lcurve(0.10)

            _, rf1, rc1, _ = zip(*points1)
            rcurve = scipy.interpolate.interp1d(rf1, rc1)
            r110 = rcurve(0.10)

            _, rf2, rc2, _ = zip(*points2)
            rcurve = scipy.interpolate.interp1d(rf2, rc2)
            r210 = rcurve(0.10)

            score1 = numpy.abs(r110/l10 - 0.90)
            score2 = numpy.abs(r210/l10 - 0.90)

            print('Resolving Rayleigh: %d %f - %d %f' % (offset1, score1,
                                                         offset2, score2))
            if score1 < score2:
                rayleighpoints = list(points1)
            else:
                rayleighpoints = list(points2)
        
    elif rayleighbounds is None:
        if not lovebounds is None:
            #
            # Special case where trough is almost halfway between reference troughs
            #
            offset1, score1, offset2, score2 = lovebounds
            print('Love between troughs:', lovedoffset, offset1, score1, offset2, score2)

            if offset1 == 0:
                points1 = list(lovepoints)
            else:
                points1 = pick(j0zeros, j1zeros, freq, lovesignal, distkm, lovephaseref,
                               args.freq_min, args.freq_max, lovedoffset + offset1)

            if offset2 == 0:
                points2 = list(lovepoints)
            else:
                points2 = pick(j0zeros, j1zeros, freq, lovesignal, distkm, lovephaseref,
                               args.freq_min, args.freq_max, lovedoffset + offset2)

            _, rf, rc, _ = zip(*rayleighpoints)
            rcurve = scipy.interpolate.interp1d(rf, rc)
            r10 = rcurve(0.10)

            _, lf1, lc1, _ = zip(*points1)
            lcurve = scipy.interpolate.interp1d(lf1, lc1)
            l110 = lcurve(0.10)

            _, lf2, lc2, _ = zip(*points2)
            lcurve = scipy.interpolate.interp1d(lf2, lc2)
            l210 = lcurve(0.10)

            score1 = numpy.abs(r10/l110 - 0.90)
            score2 = numpy.abs(r10/l210 - 0.90)

            print('Resolving Love: %d %f - %d %f' % (offset1, score1,
                                                     offset2, score2))
            if score1 < score2:
                lovepoints = list(points1)
            else:
                lovepoints = list(points2)

    else:
        print('Ambiguous/undecided')


        #
        # Special case where trough is almost halfway between reference troughs
        #
        offset1, score1, offset2, score2 = lovebounds
        print('Love between troughs:', lovedoffset, offset1, score1, offset2, score2)
        
        if offset1 == 0:
            points1 = list(lovepoints)
        else:
            points1 = pick(j0zeros, j1zeros, freq, lovesignal, distkm, lovephaseref,
                           args.freq_min, args.freq_max, lovedoffset + offset1)

        if offset2 == 0:
            points2 = list(lovepoints)
        else:
            points2 = pick(j0zeros, j1zeros, freq, lovesignal, distkm, lovephaseref,
                           args.freq_min, args.freq_max, lovedoffset + offset2)

        r10 = 1.0e9
        rok = False
        if len(rayleighpoints) >= 3 and rayleighpoints[0][1] <= 0.10:
            _, rf, rc, _ = zip(*rayleighpoints)
            rcurve = scipy.interpolate.interp1d(rf, rc)
            r10 = rcurve(0.10)
            rok = True
        
        l110 = 1.0e-3
        l1ok = False
        if len(points1) >= 3 and points1[0][1] <= 0.10:
            _, lf1, lc1, _ = zip(*points1)
            lcurve = scipy.interpolate.interp1d(lf1, lc1)
            l110 = lcurve(0.10)
            l1ok = True

        l210 = 1.0e-3
        l2ok = False
        if len(points2) >= 3 and points2[0][1] <= 0.10:
            _, lf2, lc2, _ = zip(*points2)
            lcurve = scipy.interpolate.interp1d(lf2, lc2)
            l210 = lcurve(0.10)
            l2ok = True
        
        score1 = numpy.abs(r10/l110 - 0.80)
        score2 = numpy.abs(r10/l210 - 0.80)

        if rok and l1ok and l2ok:
            print('Resolving Love: %d %f - %d %f' % (offset1, score1,
                                                     offset2, score2))
            
            if score1 < score2:
                lovepoints = list(points1)
            else:
                lovepoints = list(points2)


        offset1, score1, offset2, score2 = rayleighbounds
        print('Rayleigh between troughs:', rayleighdoffset, offset1, score1, offset2, score2)
            
        if offset1 == 0:
            points1 = list(rayleighpoints)
        else:
            points1 = pick(j0zeros, j1zeros, freq, rayleighsignal, distkm, rayleighphaseref,
                           args.freq_min, args.freq_max, rayleighdoffset + offset1)
            
        if offset2 == 0:
            points2 = list(rayleighpoints)
        else:
            points2 = pick(j0zeros, j1zeros, freq, rayleighsignal, distkm, rayleighphaseref,
                           args.freq_min, args.freq_max, rayleighdoffset + offset2)

        l10 = 1.0e-3
        lok = False
        if len(lovepoints) >= 3 and lovepoints[0][1] <= 0.10:
            _, lf, lc, _ = zip(*lovepoints)
            lcurve = scipy.interpolate.interp1d(lf, lc)
            l10 = lcurve(0.10)
            lok = True

        r110 = 1.0e9
        r1ok = False
        if len(points1) >= 3 and points1[0][1] <= 0.10:
            _, rf1, rc1, _ = zip(*points1)
            rcurve = scipy.interpolate.interp1d(rf1, rc1)
            r110 = rcurve(0.10)
            r1ok = True

        r210 = 1.0e9
        r2ok = False
        if len(points2) >= 3 and points2[0][1] <= 0.10:
            _, rf2, rc2, _ = zip(*points2)
            rcurve = scipy.interpolate.interp1d(rf2, rc2)
            r210 = rcurve(0.10)
            r2ok = True

        score1 = numpy.abs(r110/l10 - 0.90)
        score2 = numpy.abs(r210/l10 - 0.90)

        if r1ok and r2ok and lok:
            print('A: Resolving Rayleigh: %d %f - %d %f' % (offset1, score1,
                                                            offset2, score2))
            if score1 < score2:
                rayleighpoints = list(points1)
            else:
                rayleighpoints = list(points2)
            
        
    if plotting:
        peaks = []
        troughs = []
        zeros = []
        
        for (s, f, c, o) in lovepoints:
            if s > 0.0:
                peaks.append(f)
            elif s < 0.0:
                troughs.append(f)
            else:
                zeros.append(f)

                
        a = numpy.max(numpy.abs(lovesignal))
        ax[0].scatter(peaks, [a] * len(peaks), color = 'red')
        ax[0].scatter(troughs, [-a] * len(troughs), color = 'blue')
        ax[0].scatter(zeros, [0.0] * len(zeros), color = 'black')

        k = 2.0*numpy.pi*freq[indices]/lovephaseref(freq[indices])
        
        b = scipy.special.j0(k * distkm)
        mb = 0.5*numpy.max(numpy.abs(b))
        b = a/mb * b

        ax[0].plot(freq[indices], b, 'k-', linewidth = 0.5, alpha = 0.5)
        ax[0].set_ylim(-a*1.5, a*1.5)

        peaks = []
        troughs = []
        zeros = []
    
        for (s, f, c, o) in rayleighpoints:
            if s > 0.0:
                peaks.append(f)
            elif s < 0.0:
                troughs.append(f)
            else:
                zeros.append(f)

        a = numpy.max(numpy.abs(rayleighsignal))
        ax[1].scatter(peaks, [a] * len(peaks), color = 'red')
        ax[1].scatter(troughs, [-a] * len(troughs), color = 'blue')
        ax[1].scatter(zeros, [0.0] * len(zeros), color = 'black')

        k = 2.0*numpy.pi*freq[indices]/rayleighphaseref(freq[indices])
        
        b = scipy.special.j0(k * distkm)
        mb = 0.5*numpy.max(numpy.abs(b))
        b = a/mb * b

        ax[1].plot(freq[indices], b, 'k-', linewidth = 0.5, alpha = 0.5)

    if plotting:

        bx.plot(lovephaseref.x, lovephaseref.y, 'r:')
        s, f, c, o = zip(*lovepoints)
        bx.plot(f, c, 'r-')

        l = numpy.polyfit(f, c, 6)
        l = numpy.poly1d(l)
        
        ep = l(freq)
        ep = scipy.ndimage.gaussian_filter1d(ep, 3.0)
        bx.plot(freq, ep, 'k--')

        
        bx.plot(rayleighphaseref.x, rayleighphaseref.y, 'g:')
        s, f, c, o = zip(*rayleighpoints)
        bx.plot(f, c, 'g-')

        l = numpy.polyfit(f, c, 6)
        l = numpy.poly1d(l)
        
        ep = l(freq)
        ep = scipy.ndimage.gaussian_filter1d(ep, 3.0)
        bx.plot(freq, ep, 'b--')

        if args.offset != 0:

            newf = []
            newc = []
            for s, f, c, o in lovepoints:
                newo = o + args.offset
                if s == 0:
                    c = 2.0*numpy.pi*f * distkm/j0zeros[newo]
                else:
                    c = 2.0*numpy.pi*f * distkm/j1zeros[newo]

                newf.append(f)
                newc.append(c)

            bx.plot(newf, newc, 'r--')
            
            newf = []
            newc = []
            for s, f, c, o in rayleighpoints:
                newo = o + args.offset
                if s == 0:
                    c = 2.0*numpy.pi*f * distkm/j0zeros[newo]
                else:
                    c = 2.0*numpy.pi*f * distkm/j1zeros[newo]

                newf.append(f)
                newc.append(c)

            bx.plot(newf, newc, 'g--')

            

        P.show()

    if not args.output is None:
        
        f = open('%s.love' % args.output, 'w')
        for s, fr, c, o in lovepoints:

            e = estimate_error(j0zeros, j1zeros, s, fr, c, o, distkm)
            f.write('%15.9f %15.9f %d %4d %15.9f\n' % (fr, c, s, o, e))

        f.close()

        f = open('%s.rayleigh' % args.output, 'w')
        for s, fr, c, o in rayleighpoints:
            e = estimate_error(j0zeros, j1zeros, s, fr, c, o, distkm)
            f.write('%15.9f %15.9f %d %4d %15.9f\n' % (fr, c, s, o, e))

        f.close()

        
        
        
