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

def test_func(x, a, b):
    return a * np.sin(b * x) 

def plotLongitudinalMeanField(aBl):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
    fig, ax = plt.subplots(figsize=(5,4))

    ax.set_xlabel('Rotational phase')
    ax.set_ylabel('Longitudinal field (G)')

    ax.grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

    fig.suptitle("Mean longitudinal magnetic field of α2 CVn\nJuly 2022 - G. Bertrand",fontsize=9, fontweight=0, color='black' )

    for phase, (bl,err) in aBl.items():
        plt.errorbar(phase, bl, yerr=300, color="k",fmt="o", elinewidth=0.6)

    x_data = np.linspace(-0.3, 1.3, num=1000)
    #params, params_covariance = optimize.curve_fit(test_func, list(aBl.keys()), list(aBl.values()))
   
    plt.plot(x_data, test_func(x_data-0.21, 800, 6.3)-200,label='Fitted function', color="k", lw=0.6)
    plt.tight_layout(pad=1, w_pad=0, h_pad=0)
    ax.set_xlim((-0.3,1.3))
    plt.savefig('sandbox/alpha2CVn-magnetic-field-Bl.png', dpi=300)
    plt.show()

def initPlot(title):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,6))

    ax1.set_xlabel('Wavelength in Å')
    ax1.set_ylabel('I/Ic')

    ax2.set_xlabel('Wavelength in Å')
    ax2.set_ylabel('V/Ic')

    ax3.set_xlabel('Wavelength in Å')
    ax3.set_ylabel('N/Ic')

    fig.suptitle(title+"\nSkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=9, fontweight=0, color='black' )

    return [ax1, ax2, ax3]

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

def getBz(i, v, n, lambda_min, lambda_max, g= 1.0, lambda0 = 6562.8):
    v= v[(lambda_min)*u.AA:(lambda_max)*u.AA]
    n= n[(lambda_min)*u.AA:(lambda_max)*u.AA]
    i= i[(lambda_min)*u.AA:(lambda_max)*u.AA]

    # fig, ax = plt.subplots(figsize=(6,6))
    # plt.plot(v.spectral_axis, v.flux*10)
    # plt.plot(i.spectral_axis, i.flux)
    # plt.show()

    l = i.spectral_axis / 10 

    rv  = list(map(lambda x : getRv(x), i.spectral_axis))
    cc = 299792

    coeff = (-2.14 * 10 ** 11) / ((lambda0/10) * cc * g)
    bSI = 1 - i.flux  
    bSV = rv * v.flux * -1
    bSN = rv * n.flux * -1
    iSI = float(simps(bSI, rv))
    iSV = float(simps(bSV, rv))
    iSN = float(simps(bSN, rv))
    Bl = float(coeff * (iSV / iSI))
    Blerr = float(coeff * (iSN / iSI))

    return (Bl, abs(Blerr))

def plotLeftRightIVZoom(data):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'
                    
    for phase, (v, i, n, v_i, n_i, l, r) in data.items():
        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlim((6558.0,6568.5))
        #ax.set_ylim((.20,.45))

        # left
        plt.plot(l.spectral_axis, l.flux, "r-", lw=0.6, label="L = left polarization")
        # right
        plt.plot(r.spectral_axis, r.flux, "b-", lw=0.6, label="R = right polarization")
        # v 
        plt.plot(r.spectral_axis, v.flux*2+.9, "k-", lw=0.6, label="V/Ic x 2 = (L-R)/(L+R) x 2")
        # n 
        plt.plot(r.spectral_axis, n.flux*2+.9, "k--", lw=0.6, alpha=0.6, label="Null polarization spectrum")

        # phase
        ax.text(6558.35, 1.04,r'$\phi$' + ' = '+ "%.3f"%phase, size='medium')

        ax.grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

        plt.legend(loc="lower right") 

        plt.tight_layout(pad=1, w_pad=0.8, h_pad=0)
        plt.savefig('sandbox/alpha2CVn-magnetic-field-detection-individual.png', dpi=300)

        plt.show()
        break

def getStokesParam(l1, l2, r1, r2, lambda_minmaxstep):
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

    v_i = v/i
    n_i = n/i

    
    return (v, i, n, v_i, n_i, (l1+l2)/2, (r1+r2)/2)

def getPhase(jd, jd0, period):
    return round((jd-jd0) % period / period,3)

def getRv(deltalambda, lambda0 = 6562.8):
    c = 299792.458
    deltalambda = deltalambda.value
    return  (c  * (deltalambda/lambda0)) - c

#

if __name__ == '__main__':

    lambda_minmaxstep = (6546, 6582, 0.0312)

    # Eph alpha2 CVn : Farnsworth, G. 1932, ApJ, 75, 364
    phase_eph_jd0 = 2419869.720  
    phase_eph_P = 5.46939

    # Path to fits 
    path = '/Volumes/Samsung_T5/Astro/starex/corcaroli/'
    path = 'sandbox/corcaroli/'
    observations_count = 10

    # run 
    results = {}
    ax = initPlot(r"$\bf{α2}$ "+r"$\bf{CVn}$ "+"- July 2022 - Detection of mean longitudinal magnetic field - La Montagne (FR) - G. Bertrand")
    for c in range(1, observations_count+1, 1):
        l1_name = '%s-d1-h.fits'% c
        l2_name = '%s-d2-h.fits'% c
        r1_name = '%s-g1-h.fits'% c
        r2_name = '%s-g2-h.fits'% c

        l1, jd0 = createSpectrum1D(os.path.join(path, l1_name))
        l2, jd = createSpectrum1D(os.path.join(path, l2_name))
        r1, jd = createSpectrum1D(os.path.join(path, r1_name))
        r2, jd = createSpectrum1D(os.path.join(path, r2_name))

        phase = getPhase(jd0, phase_eph_jd0, phase_eph_P)
        print(phase, jd0)
        results[phase] = getStokesParam(l1,l2,r1,r2,lambda_minmaxstep)
    
    results = collections.OrderedDict(sorted(results.items(), reverse=True))
    c = 0
    stepI = 0.3
    stepV = 0.11
    for phase, (v, i, n, v_i, n_i, l, r) in results.items():
        plotSpectrum1D(ax[0], i/4 + (c*stepI))
        plotSpectrum1D(ax[1], v_i + (c*stepV))
        plotSpectrum1D(ax[2], n_i + (c*stepV))

        ax[0].text(6575, 1.06+(c*stepI),r'$\phi$' + ' = '+ "%.3f"%phase, size='medium')
        ax[1].text(6575, 0.025+(c*stepV),r'$\phi$' + ' = '+ "%.3f"%phase, size='medium')
        ax[2].text(6575, 0.025+(c*stepV),r'$\phi$' + ' = '+ "%.3f"%phase, size='medium')
        c+=1

    for i in range(0,3,1):
        ax[i].set_xlim((lambda_minmaxstep[0],lambda_minmaxstep[1]))
        ax[i].grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

    ax[0].set_ylim((0.15,4.12))
    ax[1].set_ylim((-0.12,1.12))
    ax[2].set_ylim((-0.12,1.12))

    plt.tight_layout(pad=1, w_pad=0.8, h_pad=0)
    plt.savefig('sandbox/alpha2CVn-magnetic-field-detection-%s.png' % observations_count, dpi=300)
    plt.show()

    # Individal observation
    plotLeftRightIVZoom(results)

    # longitudinal magnetic field
    # aBl = {}
    # lamb_0 = 6562.8
    # lambda_min = lamb_0 - 1.2
    # lambda_max = lamb_0 + 1.2

    # for phase, (v, i, n, v_i, n_i, l, r) in results.items():
    #     v = convolution_smooth(v, Box1DKernel(width=38))
    #     i = convolution_smooth(i, Box1DKernel(width=38))
    #     aBl[phase] = getBz(i, v, n, lambda_min, lambda_max, 1., lamb_0)

    # nBl = {}
    # for i in range(0, 4, 1):
    #     nBl[(list(aBl.keys())[i])-1] = list(aBl.values())[i]
    # for i in range(0, -4, -1):
    #     nBl[(list(aBl.keys())[i])+1] = list(aBl.values())[i]
    # aBl.update(nBl)

    # plotLongitudinalMeanField(aBl) 

    # today = date.today()
    # time = '%s-%s-%sT%s:%s:00' % ('2022', '07', '16', '22', '00')
    # t = Time(time, format='isot', scale='utc')

    # phase = getPhase(float(t.jd), phase_eph_jd0, phase_eph_P)
    # print(phase)
    # today = date.today()
    
    # exit() 