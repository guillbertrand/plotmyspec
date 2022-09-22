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
import matplotlib.pyplot as plt
from specutils import Spectrum1D,SpectralRegion
from specutils.manipulation import  LinearInterpolatedResampler, convolution_smooth
from astropy.convolution import Box1DKernel

def initPlot(title):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4,7))

    ax1.set_xlabel('Wavelength in Å')
    ax1.set_ylabel('I/Ic')

    ax2.set_xlabel('Wavelength in Å')
    ax2.set_ylabel('V/Ic')

    fig.suptitle(title+"\nSkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=9, fontweight=0, color='black' )

    return [ax1, ax2]

def plotSpectrum1D(ax, spec, color="k-"):
    return ax.plot(spec.spectral_axis * u.AA, spec.flux, color, alpha=1, lw=0.7) 

def createSpectrum1D(path):
    file = fits.open(path)  
    specdata = file[0].data
    header = file[0].header
    wcs_data = fitswcs.WCS(header={'CDELT1': round(header['CDELT1'],4), 'CRVAL1': header['CRVAL1'],
                                'CUNIT1': 'Angstrom', 'CTYPE1': 'WAVE',
                                'CRPIX1': header['CRPIX1']})
    flux= specdata * u.Jy
    s = Spectrum1D(flux=flux,  wcs=wcs_data)
    return (s, float(header['JD-OBS']))

def getStrokesParam(l1, l2, r1, r2, lambda_minmaxstep):
    lambdaMin = 6550
    lambdaMax = 6575

    # continuum normalization
    s_fit = fit_generic_continuum(l1, exclude_regions=[SpectralRegion(lambdaMin * u.AA, lambdaMax * u.AA)])
    l1   = l1 / s_fit(l1.spectral_axis)
    s_fit = fit_generic_continuum(l2, exclude_regions=[SpectralRegion(lambdaMin * u.AA, lambdaMax * u.AA)])
    l2   = l2 / s_fit(l2.spectral_axis)
    s_fit = fit_generic_continuum(r1, exclude_regions=[SpectralRegion(lambdaMin * u.AA, lambdaMax * u.AA)])
    r1   = r1 / s_fit(r1.spectral_axis)
    s_fit = fit_generic_continuum(r2, exclude_regions=[SpectralRegion(lambdaMin * u.AA, lambdaMax * u.AA)])
    r2   = r2 / s_fit(r2.spectral_axis)

    # resample spectra
    new_spectral_axis = np.arange(lambda_minmaxstep[0], lambda_minmaxstep[1], lambda_minmaxstep[2]) * u.AA
    resampler = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')
    l1 = resampler(l1, new_spectral_axis)
    l2 = resampler(l2, new_spectral_axis)
    r1 = resampler(r1, new_spectral_axis)
    r2 = resampler(r2, new_spectral_axis)

    # compute stokes IV parameters
    v = (l1 - r1) + (l2 - r2)
    i = l1 + r1 + l2 + r2
    n = (l1 - l2) + (r1 - r2)

    box1d_kernel = Box1DKernel(width=20)
    v = convolution_smooth(v, box1d_kernel)
    i = convolution_smooth(i, box1d_kernel)
    n = convolution_smooth(n, box1d_kernel)

    v_i = v/i
    n_i = n/i

    
    return (v, i, n, v_i, n_i)

def getPhase(jd, jd0, period):
    return round((jd-jd0) % period / period,3)

#

if __name__ == '__main__':

    lambda_minmaxstep = (6546, 6582,0.0312)

    # Eph alpha2 CVn : Farnsworth, G. 1932, ApJ, 75, 364
    phase_eph_jd0 = 2419869.720  
    phase_eph_P = 5.46939

    # Path to fits 
    path = '/Volumes/Samsung_T5/Astro/starex/corcaroli/'
    path = 'sandbox/corcaroli/'
    observations_count = 8

    # run 
    results = {}
    ax = initPlot(r"$\bf{α2}$ "+r"$\bf{CVn}$ "+"- July 2022 - Detection of mean longitudinal magnetic field - La Montagne (FR) - G. Bertrand")
    for c in range(1, observations_count+1, 1):
        l1_name = '%s-g1-h.fits'% c
        l2_name = '%s-g2-h.fits'% c
        r1_name = '%s-d1-h.fits'% c
        r2_name = '%s-d2-h.fits'% c

        l1, jd = createSpectrum1D(os.path.join(path, l1_name))
        l2, jd = createSpectrum1D(os.path.join(path, l2_name))
        r1, jd = createSpectrum1D(os.path.join(path, r1_name))
        r2, jd = createSpectrum1D(os.path.join(path, r2_name))

        phase = getPhase(jd, phase_eph_jd0, phase_eph_P)

        results[phase] = getStrokesParam(l1,l2,r1,r2,lambda_minmaxstep)
    
    results = collections.OrderedDict(sorted(results.items(), reverse=True))
    c = 0
    stepI = 0.3
    stepV = 0.11
    for phase, (v, i, n, v_i, n_i) in results.items():
        plotSpectrum1D(ax[0], i/4 + (c*stepI))
        plotSpectrum1D(ax[1], v_i * -1 + (c*stepV))

        ax[0].text(6576.5, 1.06+(c*stepI),r'$\phi$' + ' '+ str(phase), size='medium')
        ax[1].text(6576.5, 0.025+(c*stepV),r'$\phi$' + ' '+ str(phase), size='medium')
        c+=1

    for i in range(0,2,1):
        ax[i].set_xlim((lambda_minmaxstep[0],lambda_minmaxstep[1]))
        ax[i].grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')



    ax[0].set_xlim((6561.6,6564.2))
    ax[1].set_xlim((6561.6,6564.2))

    ax[0].set_ylim((0.15,3.7))
    ax[1].set_ylim((-0.07,1.1))

    plt.tight_layout(pad=1, w_pad=1, h_pad=0)
    plt.savefig('sandbox/alpha2CVn-magnetic-field-detection-%s.png' % observations_count, dpi=300)
    plt.show()

    
    today = date.today()
    time = '%s-%s-%sT%s:%s:00' % ('2022', '07', '13', '21', '00')
    t = Time(time, format='isot', scale='utc')

    phase = getPhase(float(t.jd), phase_eph_jd0, phase_eph_P)
    print(phase)
    today = date.today()
    
    exit() 