
import glob
import numpy
import matplotlib.pyplot as P

def safe_relative_error(v, ref):

    epsilon = 1.0e-9
    
    indices = numpy.where(numpy.abs(ref) > epsilon)[0]

    rel = numpy.zeros(v.shape)#v - ref
    rel[indices] = (v[indices] - ref[indices])/ref[indices]

    return rel

if __name__ == '__main__':

    stationpair = 'HOT05_HOT25'
    lr = 'rayleigh'
    
    files = glob.glob('Final_%s_*/opt.pred-%s' % (stationpair, lr))

    figA, ax = P.subplots()
    figB, bx = P.subplots()

    zerofilea = 'Final_%s_0/opt.pred-rayleigh' % (stationpair)
    zerofileb = 'Final_%s_0/opt.pred-love' % (stationpair)
    refa = numpy.loadtxt(zerofilea)
    refb = numpy.loadtxt(zerofileb)
    for fn in files:
        if fn == zerofilea:
            continue
        
        preda = numpy.loadtxt(fn)

        label = fn.split('/')[0][len('Final_%s_' % stationpair):]

        lfn = fn.replace('rayleigh', 'love')
        predb = numpy.loadtxt(lfn)

        rela = safe_relative_error(preda[:,5], refa[:,5])
        relb = safe_relative_error(predb[:,5], refb[:,5])

        print label, numpy.max(numpy.abs(rela)), numpy.max(numpy.abs(relb))
        
        ax.plot(preda[:,0], rela, label = label)
        bx.plot(predb[:,0], relb, label = label)

    ax.legend()
    bx.legend()

    ax.set_ylim(-0.05, 0.05)
    bx.set_ylim(-0.05, 0.05)

    ax.axhline(-0.01, color = 'red', linestyle = 'dashed')
    ax.axhline(0.01, color = 'red', linestyle = 'dashed')
    
    bx.axhline(-0.01, color = 'red', linestyle = 'dashed')
    bx.axhline(0.01, color = 'red', linestyle = 'dashed')

    P.show()
        
