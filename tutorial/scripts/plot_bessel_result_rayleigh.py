import os
import argparse

import numpy
import matplotlib.pyplot as P

from scipy.special import jacobi

def loaddispersion(fname):

    f = open(fname, 'r')
    lines = f.readlines()

    slon, slat, dlon, dlat, distkm = map(float, lines[0].split())
    freq, count, acsn, csn, N = map(float, lines[1].split())

    f, r, i, ncfr, ncfi = zip(*map(lambda x: map(float, x.split()), lines[2:]))

    spec = numpy.array(r) + numpy.array(i)*1.0j
    ncf = numpy.array(ncfr) + numpy.array(ncfi)*1.0j
    return (slon, slat, dlon, dlat, distkm, int(count)), numpy.array(f), freq, acsn, csn, spec, ncf

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', '--fits', type = str, required = True, help = 'Fits base path')

    parser.add_argument('-d', '--data', type = str, default = '../example_data', help = 'Data base path')

    parser.add_argument('--width', type = float, default = 8.0, help = 'Figure width')
    parser.add_argument('--height', type = float, default = 3.0, help = 'Figure height')
    
    args = parser.parse_args()

    if os.access(os.path.join(args.fits, 'opt.pred'), os.R_OK):
        rayleighpred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred'))
    elif os.access(os.path.join(args.fits, 'opt.pred-rayleigh'), os.R_OK):
        rayleighpred = numpy.loadtxt(os.path.join(args.fits, 'opt.pred-rayleigh'))
    else:
        raise Exception('No predictions file %s found' % os.path.join(args.fits, 'opt.pred'))

    stationpair = '_'.join(os.path.basename(args.fits.rstrip('/')).split('_')[1:3])

    rayleighdata = os.path.join(args.data, 'RayleighResponse/dispersion_%s.txt' % stationpair)

    (_, _, _, _, distkm, _), f, sample_rate, rayleigh_acsn, rayleigh_csn, rayleigh_spec, rayleigh_ncf = loaddispersion(rayleighdata)

    figB, bx = P.subplots()
    figB.set_size_inches((args.width, args.height))
    figB.set_tight_layout(True)

    #
    # Modulated Bessel is column 5, Raw Bessel is column 3, Envelope 4
    #
    colindex = 5
    #colindex = 3

    indices = numpy.where(rayleighpred[:,1] > 0.0)[0]

    bx.set_title('Rayleigh')
    bx.plot(rayleighpred[indices,0], rayleighpred[indices,colindex], 'r-', linewidth = 1, zorder = 100)
    bx.plot(f, numpy.real(rayleigh_ncf), linestyle = 'solid', color = 'grey', linewidth = 2, zorder = 50)

    bx.set_xlim(0, 0.4)

    bx.set_xlabel('Frequency (Hz)')

    P.show()
