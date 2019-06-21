
import os

import argparse

import numpy
import matplotlib.pyplot as P
import matplotlib.gridspec

import scipy.signal
from scipy.special import jacobi

def autosigma(distkm):
    return 0.4418954398283702/numpy.sqrt(distkm) - 0.007296666006375768

def loaddispersion(fname):

    f = open(fname, 'r')
    lines = f.readlines()

    slon, slat, dlon, dlat, distkm = map(float, lines[0].split())
    freq, count, acsn, csn, N = map(float, lines[1].split())

    f, r, i, ncfr, ncfi = zip(*map(lambda x: map(float, x.split()), lines[2:]))

    spec = numpy.array(r) + numpy.array(i)*1.0j
    ncf = numpy.array(ncfr) + numpy.array(ncfi)*1.0j
    return (slon, slat, dlon, dlat, distkm, int(count)), numpy.array(f), freq, acsn, csn, spec, ncf

def mkftan(f, sample_rate, spec, period_min, period_max, vmin, vmax, vsample, distkm, sigma, period = False):

    vaxis = numpy.linspace(vmin, vmax, vsample)
    vtaxis = distkm/vaxis
    N = len(spec) - 1
    t = (numpy.arange(N, dtype = 'float') + 0.5)/sample_rate

    causal_vimg = numpy.zeros((vsample, N + 1))
    acausal_vimg = numpy.zeros((vsample, N + 1))
    
    for i in range(1, len(f)):

        c = f[i]
        period = 1.0/f[i]
        if (period >= period_min and period <= period_max):
            g = numpy.exp(-(f - c)**2/(2.0 * sigma**2))

            fspec = g*spec

            pt = numpy.fft.irfft(fspec)
            ft = numpy.abs(scipy.signal.hilbert(pt))
            acausal = ft[::-1][:N]
            causal = ft[:N]

            vinterp = scipy.interpolate.interp1d(t[:N], causal)
            causal_vimg[:,i] = vinterp(vtaxis)
            vinterp = scipy.interpolate.interp1d(t[:N], acausal)
            acausal_vimg[:,i] = vinterp(vtaxis)

    return vtaxis, causal_vimg, acausal_vimg


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fits', type = str, required = True, help = 'Fits base path')

    parser.add_argument('-d', '--data', type = str, default = '../example_data', help = 'Data base path')

    parser.add_argument('-S', '--sigma', type = float, default = 0.02, help = 'Gaussian filter sigma')
    parser.add_argument('-A', '--auto-sigma', action = 'store_true', default = False, help = 'Automatic sigma')

    parser.add_argument('-v', '--vmin', type = float, default = 1.0, help = 'Minimum velocity')
    parser.add_argument('-V', '--vmax', type = float, default = 5.0, help = 'Maximum velocity')
    parser.add_argument('-N', '--vsample', type = int, default = 500, help = 'Velocity FTAN samples')

    parser.add_argument('-p', '--period-min', type = float, default = 2.0, help = 'Min period')
    parser.add_argument('-P', '--period-max', type = float, default = 40.0, help = 'Max period')

    parser.add_argument('--pdf', type = str, default = None, help = 'PDF')

    parser.add_argument('--cmap', type = str, default = 'Greys', help = 'Colour map')
    
    parser.add_argument('--width', type = float, default = 8.0, help = 'Figure width')
    parser.add_argument('--height', type = float, default = 3.0, help = 'Figure height')
    
    args = parser.parse_args()

    lovepred = numpy.loadtxt('%s/opt.pred-love' % args.fits)
    rayleighpred = numpy.loadtxt('%s/opt.pred-rayleigh' % args.fits)

    stationpair = '_'.join(os.path.basename(args.fits.rstrip('/')).split('_')[1:3])

    lovedata = os.path.join(args.data, 'LoveResponse/dispersion_%s.txt' % stationpair)
    rayleighdata = os.path.join(args.data, 'RayleighResponse/dispersion_%s.txt' % stationpair)

    (_, _, _, _, distkm, _), f, sample_rate, love_acsn, love_csn, love_spec, love_ncf = loaddispersion(lovedata)
    
    (_, _, _, _, distkm, _), f, sample_rate, rayleigh_acsn, rayleigh_csn, rayleigh_spec, rayleigh_ncf = loaddispersion(rayleighdata)

    print(distkm, args.sigma, distkm/args.sigma, args.sigma/distkm)

    sigma = autosigma(distkm)
    print('Auto sigma:', sigma)

    vtaxis, causal_love, acausal_love = mkftan(f, sample_rate, love_spec,
                                               args.period_min, args.period_max,
                                               args.vmin, args.vmax, args.vsample, distkm, sigma)

    vtaxis, causal_rayleigh, acausal_rayleigh = mkftan(f, sample_rate, rayleigh_spec,
                                                       args.period_min, args.period_max,
                                                       args.vmin, args.vmax, args.vsample, distkm, sigma)



    figA, ax = P.subplots()
    figA.set_tight_layout(True)
    figA.set_size_inches((args.width, args.height))

    mlove = acausal_love
    thresh = 0.1 * numpy.max(mlove)
        
    r, c = mlove.shape
            
            
    ax.contourf(mlove, extent = [f[0], f[-1], args.vmin, args.vmax],
                cmap = args.cmap)
    indices = numpy.where(lovepred[:,1] > 0.0)[0]
    ax.plot(lovepred[indices,0], lovepred[indices,1]/1.0e3, 'k-')

    ax.set_xlim(0.0, 0.4)
    ax.set_ylim(args.vmin, args.vmax)
    ax.set_ylabel('Group velocity (km/s)')
    ax.set_xlabel('Frequency (Hz)')

    figB, bx = P.subplots()
    figB.set_tight_layout(True)
    figB.set_size_inches((args.width, args.height))

    mrayleigh = acausal_rayleigh
        
    thresh = 0.1 * numpy.max(mrayleigh)
        
    r, c = mrayleigh.shape
            
            
    bx.contourf(mrayleigh, extent = [f[0], f[-1], args.vmin, args.vmax],
                cmap = args.cmap)
    indices = numpy.where(rayleighpred[:,1] > 0.0)[0]
    bx.plot(rayleighpred[indices,0], rayleighpred[indices,1]/1.0e3, 'k-')

    bx.set_xlim(0.0, 0.4)
    bx.set_ylim(args.vmin, args.vmax)
    bx.set_ylabel('Group velocity (km/s)')
    bx.set_xlabel('Frequency (Hz)')

    if args.pdf is None:
        P.show()
    else:
        figA.savefig('%s_love.pdf' % args.pdf, format = 'PDF')
        figB.savefig('%s_rayleigh.pdf' % args.pdf, format = 'PDF')
    
    




