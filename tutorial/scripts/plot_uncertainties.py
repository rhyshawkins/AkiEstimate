
import matplotlib.pyplot as P

import os
import numpy
import scipy.interpolate
import scipy.ndimage.filters

import scipy.stats

import argparse

def uncertainty(path):

    #
    # Posterior uncertainty on model CM = (G^T C_d^-1 G + C_m^-1)^-1
    #
    # Posterior uncertainty on phase = Jc CM Jc.T
    # Posterior uncertainty on group = JU CM JU.T
    #


    #
    # Cm was saved unsquared
    #
    Cm = numpy.loadtxt(os.path.join(path, 'opt.Cm'))**2
    
    Cd_love = numpy.loadtxt(os.path.join(path, 'opt.love_Cd'))
    Cd_rayleigh = numpy.loadtxt(os.path.join(path, 'opt.rayleigh_Cd'))

    G_love = numpy.loadtxt(os.path.join(path, 'opt.love_G'))
    G_rayleigh = numpy.loadtxt(os.path.join(path, 'opt.rayleigh_G'))

    Nd = Cd_love.size + Cd_rayleigh.size
    Nd_love, Nm = G_love.shape
    Nd_rayleigh, _ = G_rayleigh.shape

    G = numpy.zeros((Nd, Nm))
    G[:Nd_love, :] = G_love
    G[Nd_love:, :] = G_rayleigh
    
    Cd = numpy.zeros((Nd,))

    Cd[:Nd_love] = Cd_love
    Cd[Nd_love:] = Cd_rayleigh

    CMinv = (G.T.dot(numpy.diag(1.0/Cd)).dot(G) + 
             numpy.diag(1.0/Cm))

    CM = numpy.linalg.inv(CMinv)
    
    Jc_love = numpy.loadtxt(os.path.join(path, 'opt.love_Jc'))
    Jc_love = Jc_love[::-1,:]
    
    Jc_rayleigh = numpy.loadtxt(os.path.join(path, 'opt.rayleigh_Jc'))
    Jc_rayleigh = Jc_rayleigh[::-1,:]
                       
    JU_love = numpy.loadtxt(os.path.join(path, 'opt.love_JU'))
    JU_love = JU_love[::-1,:]
    
    JU_rayleigh = numpy.loadtxt(os.path.join(path, 'opt.rayleigh_JU'))
    JU_rayleigh = JU_rayleigh[::-1,:]

    Cc_love = Jc_love.dot(CM).dot(Jc_love.T)
    Cc_rayleigh = Jc_rayleigh.dot(CM).dot(Jc_rayleigh.T)

    CU_love = JU_love.dot(CM).dot(JU_love.T)
    CU_rayleigh = JU_rayleigh.dot(CM).dot(JU_rayleigh.T)

    return CM, Cc_love, Cc_rayleigh, CU_love, CU_rayleigh

def eig_covariance_projection(cov, dof):
    r, c = cov.shape

    F = scipy.stats.chi2.ppf(0.95, df = dof)
    
    proj = numpy.zeros((r,))
    
    w, v = numpy.linalg.eigh(cov)

    for i in range(c):

        if w[i] > 0.0:

            p = numpy.abs(v[:,i] * numpy.sqrt(w[i] * F))
            proj = numpy.fmax(proj, 2.0*p)

    return proj

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fits', type = str, required = True, help = 'Results path')

    parser.add_argument('--width', type = float, default = 8.0, help = 'Width')
    parser.add_argument('--height', type = float, default = 3.0, help = 'Height')

    parser.add_argument('--zoom', action = 'store_true', default = False, help = 'Zoom box')
    
    parser.add_argument('--pdf', type = str, default = None, help = 'PDF output')

    args = parser.parse_args()
    
    CM, Cc_love, Cc_rayleigh, CU_love, CU_rayleigh = uncertainty(args.fits)

    dofM, _ = CM.shape
    dofD, _ = Cc_rayleigh.shape

    #
    # Phase Love
    #
    figA, ax = P.subplots()
    figA.set_tight_layout(True)
    figA.set_size_inches((args.width, args.height))
    
    pred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred-love'))
    indices = numpy.where(pred[:,1] > 0.0)[0]
    
    f = pred[indices,0]
    U = pred[indices,1]
    c = pred[indices,2]

    ax.plot(f, c/1.0e3, 'k-')

    r = eig_covariance_projection(Cc_love, dofD)

    ax.fill_between(f, (c + r)/1.0e3, (c - r)/1.0e3, color = 'grey', alpha = 0.5)
    ax.set_ylim(2, 4.5)
    ax.set_xlim(0, 0.4)

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Phase velocity (km/s)')

    #
    # Phase Rayleigh
    #
    figB, bx = P.subplots()
    figB.set_tight_layout(True)
    figB.set_size_inches((args.width, args.height))
    
    pred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred-rayleigh'))
    indices = numpy.where(pred[:,1] > 0.0)[0]
    
    f = pred[indices,0]
    U = pred[indices,1]
    c = pred[indices,2]

    bx.plot(f, c/1.0e3, 'k-')

    r = eig_covariance_projection(Cc_rayleigh, dofD)

    bx.fill_between(f, (c + r)/1.0e3, (c - r)/1.0e3, color = 'grey', alpha = 0.5)
    bx.set_ylim(2, 4.5)
    bx.set_xlim(0, 0.4)

    bx.set_xlabel('Frequency (Hz)')
    bx.set_ylabel('Phase velocity (km/s)')

    #
    # Group Love
    #
    figC, cx = P.subplots()
    figC.set_tight_layout(True)
    figC.set_size_inches((args.width, args.height))
    
    pred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred-love'))
    indices = numpy.where(pred[:,1] > 0.0)[0]
    
    f = pred[indices,0]
    U = pred[indices,1]
    c = pred[indices,2]

    cx.plot(f, U/1.0e3, 'k-')

    r = eig_covariance_projection(CU_love, dofD)

    cx.fill_between(f, (U + r)/1.0e3, (U - r)/1.0e3, color = 'grey', alpha = 0.5)
    cx.set_ylim(1, 5)
    cx.set_xlim(0, 0.4)

    cx.set_xlabel('Frequency (Hz)')
    cx.set_ylabel('Group velocity (km/s)')

    #
    # Group Rayleigh
    #
    figD, dx = P.subplots()
    figD.set_tight_layout(True)
    figD.set_size_inches((args.width, args.height))
    
    pred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred-rayleigh'))
    indices = numpy.where(pred[:,1] > 0.0)[0]
    
    f = pred[indices,0]
    U = pred[indices,1]
    c = pred[indices,2]

    dx.plot(f, U/1.0e3, 'k-')

    r = eig_covariance_projection(CU_rayleigh, dofD)

    dx.fill_between(f, (U + r)/1.0e3, (U - r)/1.0e3, color = 'grey', alpha = 0.5)
    dx.set_ylim(1, 5)
    dx.set_xlim(0, 0.4)

    dx.set_xlabel('Frequency (Hz)')
    dx.set_ylabel('Group velocity (km/s)')

    if args.pdf is None:
        P.show()
    else:
        figA.savefig('%s_love_phase.pdf' % args.pdf, format = 'PDF')
        figB.savefig('%s_rayleigh_phase.pdf' % args.pdf, format = 'PDF')
        figC.savefig('%s_love_group.pdf' % args.pdf, format = 'PDF')
        figD.savefig('%s_rayleigh_group.pdf' % args.pdf, format = 'PDF')
    

    
