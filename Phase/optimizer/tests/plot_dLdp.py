
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
    files = glob.glob('Final_%s_*/opt.dLdp' % stationpair)

    figA, ax = P.subplots()
    figB, bx = P.subplots()

    zerofile = 'Final_%s_0/opt.dLdp' % stationpair
    ref = numpy.loadtxt(zerofile)
    
    for fn in files:
        if fn == zerofile:
            continue
        
        dLdp = numpy.loadtxt(fn)

        label = fn.split('/')[0][len('Final_%s_' % stationpair):]

        ax.plot(safe_relative_error(dLdp[:,0], ref[:,0]), label = label)
        bx.plot(safe_relative_error(dLdp[:,1], ref[:,1]), label = label)

    ax.legend()
    bx.legend()

    ax.set_ylim(-0.05, 0.05)
    bx.set_ylim(-0.05, 0.05)

    ax.axhline(-0.01, color = 'red', linestyle = 'dashed')
    ax.axhline(0.01, color = 'red', linestyle = 'dashed')
    
    bx.axhline(-0.01, color = 'red', linestyle = 'dashed')
    bx.axhline(0.01, color = 'red', linestyle = 'dashed')

    P.show()
        
