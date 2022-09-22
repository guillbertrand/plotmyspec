from cProfile import label
import os
import collections
import numpy as np
from astropy import units as u 
from astropy.io import fits
from astropy.time import Time
from astropy.modeling import models, fitting
from specutils.fitting import fit_generic_continuum
from datetime import date
import astropy.wcs as fitswcs 
from scipy.integrate import simps
import matplotlib.pyplot as plt
from scipy import misc, optimize
from specutils import Spectrum1D,SpectralRegion
from specutils.manipulation import  LinearInterpolatedResampler, convolution_smooth
from astropy.convolution import Box1DKernel
from scipy import optimize

def initPlot(title):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
    f = plt.figure(figsize=(8,8))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(223)
    ax3 = plt.subplot(122)

    ax1.set_xlabel('Wavelength in Å')
    ax1.set_ylabel('Relative intensity')

    ax2.set_xlabel('Wavelength in Å')
    ax2.set_ylabel('Relative intensity')

    ax3.set_xlabel('Wavelength in Å')
    ax3.set_ylabel('Relative intensity')

    plt.suptitle(title+"\n2022/08/27 - La Montagne (FR) - G. Bertrand\nSkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=9, fontweight=0, color='black' )

    return [ax1, ax2, ax3]

def plotSpectrum1D(ax, spec, label, color="k-", alpha=1, shift=0):
    if color == None:
        return ax.plot(spec.spectral_axis +shift * u.AA, spec.flux, alpha=alpha, lw=0.7, label=label) 
    else:
        return ax.plot(spec.spectral_axis +shift * u.AA, spec.flux, color, alpha=alpha, lw=0.7, label=label) 

def createSpectrum1D(path):
    file = fits.open(path)  
    specdata = file[0].data
    header = file[0].header
    print(specdata)
    wcs_data = fitswcs.WCS(header={'CDELT1': header['CDELT1'], 'CRVAL1': header['CRVAL1'],
                                'CUNIT1': 'Angstrom', 'CTYPE1': 'WAVE',
                                'CRPIX1': header['CRPIX1']})


        # Get first pixel reference
    xRef = header['CRPIX1'] - 1
    #Get length of data axis1
    xlength = header['NAXIS1']
    #Get Wavelength pixel step (Angstrom) and convert
    lambdaStep = header['CDELT1']
    lambda1 = header['CRVAL1'] - lambdaStep * xRef
    lambda2 = header['CRVAL1'] + lambdaStep * (xlength - xRef)

    wavelength = np.arange(lambda1, lambda2, lambdaStep) * u.AA
    flux= specdata * u.Jy
    s = Spectrum1D(flux=flux,  spectral_axis=wavelength)
    return (s, float(header['JD-OBS']))

def getRv(deltalambda, lambda0 = 6562.8):
    c = 299792.458
    return  (c  * (deltalambda/lambda0)) 

#

if __name__ == '__main__':

    lambda_minmaxstep = [(5875.0, 5911.0),(5875.0, 5911.0),(5876.0, 5909.0)]

    # Path to fits 
    io = 'sandbox/io/io_shifted.fits'
    europa = 'sandbox/io/europa_shifted.fits'
    ganymede = 'sandbox/io/ganymede.fits'

    io0 = 'sandbox/io/_io_20220827_946.fits'
    europa0 = 'sandbox/io/_europa_20220827_971.fits'
    ganymede0 = 'sandbox/io/_ganymede_20220827_923.fits'

    # run 
    results = {}
    
    spectrum_title = ''
    main= "High resolution spectra of Io's neutral sodium cloud"
    split_oname = main.split(' ')
    for w in split_oname:
        spectrum_title += r"$\bf{%s}$ " % (w)

    ax = initPlot(spectrum_title)

    io, jd0 = createSpectrum1D(io)
    europa, jd0 = createSpectrum1D(europa)
    ganymede, jd0 = createSpectrum1D(ganymede)
    diff_e_g, jd0 = createSpectrum1D('sandbox/io/diff_europa_ganymede.fits')
    diff_i_g, jd0 = createSpectrum1D('sandbox/io/diff_io_ganymede.fits')
    jupiter0, jd0 = createSpectrum1D('sandbox/io/jupiter_full.fits')

    io0, jd0 = createSpectrum1D(io0)
    europa0, jd0 = createSpectrum1D(europa0)
    ganymede0, jd0 = createSpectrum1D(ganymede0)

    plotSpectrum1D(ax[0], jupiter0, 'Jupiter', None)
    plotSpectrum1D(ax[0], ganymede0, 'Ganymede', None)
    plotSpectrum1D(ax[0], europa0 , 'Europa', None)
    plotSpectrum1D(ax[0], io0, 'Io', None)
    ax[0].legend() 
    ax[0].set_title("Raw spectra of Io, Ganymede and Europa",fontweight="bold", size=8)

    s=0.55
    plotSpectrum1D(ax[1], ganymede, 'Ganymede', None, shift=s)
    plotSpectrum1D(ax[1], europa , 'Europa', None,shift=s)
    plotSpectrum1D(ax[1], io, 'Io', None,shift=s)
    ax[1].legend() 
    ax[1].set_title("Same with H2O telluric lines removed\n and respective radial velocity corrected",fontweight="bold", size=8)

    plotSpectrum1D(ax[2], ganymede-0.1-.8, 'Ganymede',shift=s)
    plotSpectrum1D(ax[2], europa+0.9-.8 , 'Europa',shift=s)
    plotSpectrum1D(ax[2], io+.8+.2, 'Io',shift=s)
    plotSpectrum1D(ax[2], (diff_i_g+4.5), 'Ganymede',shift=s)
    plotSpectrum1D(ax[2], (diff_e_g+.5+.9+.9+.9+.4), 'Ganymede',shift=s)
    
    lines = [('D1', 5895.924, .5, 0.18),('D2', 5889.950, .5, 0.18)]
    for line in lines:
        n = ''
        n, lam, offset_x, offset_y = line
        ax[0].axvline(x=float(lam), color='black', linestyle='--', linewidth=0.7, alpha=0.6)
        ax[1].axvline(x=float(lam), color='black', linestyle='--', linewidth=0.7, alpha=0.6)
        ax[2].axvline(x=float(lam), color='black', linestyle='--', linewidth=0.7, alpha=0.6)
        ax[0].text(float(lam)+float(offset_x), float(offset_y), n, va='bottom', color='k',size='7')
        ax[1].text(float(lam)+float(offset_x), float(offset_y), n, va='bottom', color='k',size='7')
        ax[2].text(float(lam)+float(offset_x), float(offset_y)-1.1, n, va='bottom', color='k',size='7')

    ax[2].text(5876.5, .75-.8, 'Ganymede', size='medium')
    ax[2].text(5876.5, .75+1-.8, 'Europa', size='medium')
    ax[2].text(5876.5, .7+2-.8, 'Io', size='medium')
    ax[2].text(5876.5, 4-.7, 'Difference of Europa and Ganymede spectra', size='medium')
    ax[2].text(5876.5, .3+4, 'Difference of Io and Ganymede spectra', size='medium')

    for i in range(0,3,1):
        ax[i].set_xlim((lambda_minmaxstep[i][0],lambda_minmaxstep[i][1]))
        ax[i].grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

    #ax[0].set_ylim((0.15,4.12))
    ax[2].set_ylim((-1,6.2))

    plt.tight_layout(pad=1.5, w_pad=1.5, h_pad=1.5)
    plt.savefig('sandbox/io-sodium-cloud.png', dpi=300)
    plt.show()

    print(getRv((5890.82-5889.950), 5889.950))


