from cProfile import label
import os
import numpy as np
from astropy import units as u #units
from astropy.io import fits
from astropy.time import Time
from datetime import date
import astropy.wcs as fitswcs #wcs
import matplotlib.pyplot as plt
from specutils import Spectrum1D, SpectralRegion 
from specutils.manipulation import extract_region
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler

def initPlot(title):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel('Wavelength in Å')
    ax.set_ylabel('Relative intensity')

    #Grid configuration
    ax.grid(color='grey', alpha=0.8, linestyle='-.', linewidth=0.2, axis='both') 

    plt.suptitle(title,fontsize=9, fontweight=0, color='black' )
    plt.title("SkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=8, fontweight=0, color='black')

    return ax, ax.twinx()

def plotSpectrum1D(ax, spec, color="k-"):
    return ax.plot(spec.spectral_axis * u.AA, spec.flux, color, alpha=1, lw=0.6) 

def createSpectrum1D(path):
    #open & load spectrum file
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
    new_spectral_axis = np.arange(lambda_minmaxstep[0], lambda_minmaxstep[1], lambda_minmaxstep[2]) * u.AA
    resampler = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')
    l1 = resampler(l1, new_spectral_axis)
    l2 = resampler(l2, new_spectral_axis)
    r1 = resampler(r1, new_spectral_axis)
    r2 = resampler(r2, new_spectral_axis)
    v = (r1 - l1) + (r2 - l2)
    i = l1 + r1 + l2 + r2
    n = (l1 - l2) + (r1 - r2)
    v_i = v/i
    n_i = n/i
    return (v, i, n, v_i, n_i)

def getPhase(jd, jd0, period):
    return round((jd-jd0) % period / period,3)

#

path = '/Volumes/Samsung_T5/Astro/starex/corcaroli/'

l1_name = '2_g1.fits'
l2_name = '2_g2.fits'
r1_name = '2_d1.fits'
r2_name = '2_d2.fits'

lambda_minmaxstep = (6546, 6582,0.0312)
phase_eph_jd0 = 2419869.720  
phase_eph_P = 5.46939
''' 
today = date.today()
time = '%s-%s-%sT%s:%s:00' % ('2022', '07', '02', '21', '00')
t = Time(time, format='isot', scale='utc')

phase = getPhase(float(t.jd), phase_eph_jd0, phase_eph_P)
print(phase)
exit() '''
l1, jd = createSpectrum1D(os.path.join(path, l1_name))
l2, jd = createSpectrum1D(os.path.join(path, l2_name))
r1, jd = createSpectrum1D(os.path.join(path, r1_name))
r2, jd = createSpectrum1D(os.path.join(path, r2_name))
(v, i, n, v_i, n_i) = getStrokesParam(l1,l2,r1,r2,lambda_minmaxstep)

phase = getPhase(jd, phase_eph_jd0, phase_eph_P)

(ax, ax2) = initPlot("2022/06/30 - Detection of mean longitudinal magnetic field on α2 CVn")
l1, = plotSpectrum1D(ax, i/4 + 0.2, 'k-')
l2, = plotSpectrum1D(ax2, v_i, 'b-')
l3, = plotSpectrum1D(ax2, n_i, 'r--')
plt.tight_layout(pad=1, w_pad=0, h_pad=0)
plt.text(6579, 0.03,r'$\phi$' + ' '+ str(phase), size='medium')
ax.set_xlim((lambda_minmaxstep[0],lambda_minmaxstep[1]))
ax.set_ylim((0,1.4))
ax2.set_ylim((-0.1,0.55))
plt.legend([l1, l2, l3], ['I/Ic', 'V/Ic', 'N/Ic'])
plt.savefig('sandbox/alpha2CVn-magnetic-field-detection-2.png', dpi=300)
plt.show()

