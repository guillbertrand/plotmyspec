import logging
import re
import os
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from orbital import utilities
from astropy.io import fits
import specutils
from specutils import Spectrum1D,SpectralRegion
import astropy.wcs as fitswcs #wcs
from specutils.fitting import find_lines_derivative, fit_lines, estimate_line_parameters
from specutils.manipulation import extract_region
from specutils.analysis import centroid
from astropy.modeling import models, fitting
from specutils.fitting import fit_generic_continuum


from datetime import date

def binarySystemObservation(obj, filename, period, jd0, step=0.05):
    """
    Returns the dates of the next possible observations 
    for all phases (according to the step) of the period (P)

    Parameters:
    obj (str): Object name
    filename (str): Path of your time.lst file from ISIS
                    example time.lst : 
                            @m9 -0.318175
                            @m1 0
                            @m2 0.194897
                            @m3 0.340229
 
    period (float): Perdiod P of the object
    jd0 (float): First julian date of observation
    step: size of step for increment the phase 

    Returns:
    Log results 

    """
    logging.info('\U00002728 Start scan for %s period (%s days)' % (obj, period))
    data = []
    results = []
    with open(filename, newline='\n') as f:
        for line in f:
            p = float(line.rstrip().split(' ')[1])
            if(p>=0):
                data.append(p)
        
    for i in np.arange(0,1,step):
        res = None
        for d in data:
            if(d==i):
                res = d
                break

            ip = i + step
            im = i - step
            if(ip>1):
                ip=im
                im=(i + step)%1
            if(im<0):
                im=ip
                ip=(i - step)%1
            if (d > im and d < ip):
                res = d
                break
        
        t = findNextDateByPhase(period, i, jd0, step)
        if(res != None):
            results.append('\U00002705 Phase = %f : next date=> %s (ok : %s)' % ( round(i,2),t.iso, res))
        else:
            results.append('\U0000274C Phase = %f : next date => %s' % (round(i,2), t.iso))
    
    for r in results:
        logging.info(r)
            
def findNextDateByPhase(period, phase, jd0, step):
    data = []
    today = date.today()
    time = '%s-%s-%sT%s:%s:00' % (today.strftime("%Y"), today.strftime("%m"), today.strftime("%d"), 23, 00)
    t = Time(time, format='isot', scale='utc')
    res = None
    while not res:
        p =  (t.jd-jd0) % period / period
        if(p > (phase - step) and p < phase + step):
            return t
        t = t + 1*u.d 

def getPhase(jd0, period, jd):
    return (jd-jd0) % period / period

def getBinSysData(specs, period):
    obs = {}
    for s in specs:
        f = fits.open(s)
        header = f[0].header
        jd = header['JD-OBS']
        obs[jd] = {'fits':s, 'centroid':float(header['S_CAL'])}

    jd0 = min(obs.keys())
    for jd in obs:
        phase = getPhase(float(jd0), period, float(jd))
        obs[jd]['phase'] = phase
    return obs

def findCentroid(s):
    #open & load spectrum file
    file = fits.open(s['fits'])  
    specdata = file[0].data
    header = file[0].header

    #make WCS object
    wcs_data = fitswcs.WCS(header={'CDELT1': header['CDELT1'], 'CRVAL1': header['CRVAL1'],
                                'CUNIT1': header['CUNIT1'], 'CTYPE1': 'WAVE',
                                'CRPIX1': header['CRPIX1']})

    #set flux units
    flux= specdata * u.Jy

    #create a Spectrum1D object with specutils
    spec = Spectrum1D(flux=flux, wcs=wcs_data)

    #prepare data
    x = spec.spectral_axis
    y = spec.flux

    #fitting continuum (with exclude region between 3700A and 4000A)
    g1_fit = fit_generic_continuum(spec, exclude_regions=[SpectralRegion(6400 * u.AA, 6600 * u.AA)])
    y_continuum_fitted = g1_fit(x)
    spec_normalized = spec / y_continuum_fitted

    #Line fitting with Derivative technique
    lines = find_lines_derivative(spec_normalized, flux_threshold=-0.4)
    line = None
    for l in lines:
        if(l['line_type'] == 'absorption'):
            line = l['line_center'].value*1e10
            break


    #
    sr =  SpectralRegion((line-0.8)*u.AA, (line+0.8)*u.AA)

    # spec_normalized = extract_region(spec_normalized, sr)
    # f, ax = plt.subplots()  
    # ax.plot(spec_normalized.spectral_axis, spec_normalized.flux)   
    # plt.grid(True)  
    
    return centroid(spec_normalized, sr).value*1e10 #+0.15


def getRadialVelocity(t, t0, P, K, e, w, v0):
    # Mean anomaly
    M = 2 * np.pi *  ((t - t0)%1) / 1   
    # Eccentric anomaly
    E = utilities.eccentric_anomaly_from_mean(e,M, tolerance=0.00001)   
    # True anomaly
    f = utilities.true_anomaly_from_eccentric(e, E) 
    return (K * (e * np.cos(w) + np.cos(w + f)) + v0 ) 

def getRv(deltalambda, lambda0 = 6562.8):
    return  (299792.458 * ((deltalambda/lambda0) ** 2 - 1) / ((deltalambda/lambda0) ** 2 + 1))

def plotRadialVelocityCurveSB1(filename, P, v0, K, e, w, jd0):
    xp = []
    xjd = []
    yrv = []
    yerr = []
    with open(filename, newline='\n') as f:
        for line in f:
            p = line.rstrip().split(' ')
            xjd.append(float(p[0]))
            xp.append(float(p[1]))
            lambda0 = 6562.8
            deltalambda = (float(p[2])+float(p[3])+float(p[4]))/3
            rv = getRv(deltalambda)
            yrv.append(rv) 
            deltalambda_min = min(np.float64(p[2:]))
            rv_min = getRv(deltalambda_min)
            deltalambda_max = max(np.float64(p[2:]))
            rv_max = getRv(deltalambda_max)
            err = rv_max - rv_min
            yerr.append(err)
    
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'

    model_x = np.arange(-0.5,1.01, 0.005)
    model_y = list(map(lambda x: getRadialVelocity(x,jd0,P,K,e,w,v0), model_x))
    fig, ax =  plt.subplots(figsize=(11,6))

    #Add Graph title
    plt.suptitle(r"$\bf{HD123299}$" + " - α Dra - Phased radial-velocities - 4 observations collected from April to June 2022",fontsize=9, fontweight=0, color='black' )
    plt.title("SkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=8, fontweight=0, color='black')

    #Add X axis label
    ax.set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)

    #Add Y axis label
    ax.set_ylabel('Radial velocity [km/s]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    
    plt.errorbar(xp, yrv, yerr=yerr, color="black",fmt="o",lw=0.8)
    plt.plot(model_x, model_y, alpha=1, color="red", lw=0.5, label='"official" curve parameters - 2017 - arXiv:1707.05090')

    ax.grid(color='grey', alpha=0.4, linestyle='-', linewidth=0.5, axis='both')

    plt.legend() 

    plt.tight_layout(pad=1, w_pad=0, h_pad=0)

    plt.xticks(np.arange(-0.5, 1.1, 0.1))
    plt.yticks(np.arange(-50, 60, 10))
    plt.savefig('sandbox/alphadra/hd123299-phased-radial-velocities.png', dpi=300)
    plt.show()  

def plotRadialVelocityCurveSB1_2(specs, P, v0, K, e, w, jd0):
    xp = []
    xjd = []
    yrv = []
    yerr = []
   
    for jd, s in specs.items():
        xjd.append(jd)
        xp.append(s['phase'])
        yrv.append(s['rv']) 
    
    xp.append(1)
    xjd.append(min(xjd))
    yrv.append(specs[min(xjd)]['rv'])
    
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'

    model_x = np.arange(-0.1,1.01, 0.005)
    model_y = list(map(lambda x: getRadialVelocity(x,jd0,P,K,e,w,v0), model_x))
    fig, ax =  plt.subplots(figsize=(11,6))

    #Add Graph title
    plt.suptitle(r"$\bf{HD123299}$" + " - α Dra - Phased radial-velocities - %s observations collected from April to June 2022" % len(specs),fontsize=9, fontweight=0, color='black' )
    plt.title("SkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=8, fontweight=0, color='black')

    #Add X axis label
    ax.set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)

    #Add Y axis label
    ax.set_ylabel('Radial velocity [km/s]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    
    plt.plot(xp, yrv, "o", color="black", lw=0.8)
    plt.plot(model_x, model_y, alpha=1, color="red", lw=0.5, label='"official" curve parameters - 2017 - arXiv:1707.05090')

    ax.grid(color='grey', alpha=0.4, linestyle='-', linewidth=0.5, axis='both')

    plt.legend() 

    plt.tight_layout(pad=1, w_pad=0, h_pad=0)

    plt.xticks(np.arange(-0.1, 1.1, 0.1))
    plt.yticks(np.arange(-50, 60, 10))
    plt.savefig('sandbox/alphadra/hd123299-phased-radial-velocities-auto.png', dpi=300)
    plt.show()  
      

if __name__ == '__main__':
    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 Planning observations of a binary system - Start \U0001F680')

    # Run for mizar
    filename = 'D:\\ASTRO\\Starex\\mizar\\time2.lst'
    res = binarySystemObservation('mizar', filename, 20.53835, 2459720.381331, 0.045)

    # Run for alpha dra
    filename = 'D:\\ASTRO\\Starex\\alphadra\\time2.lst'
    res = binarySystemObservation('alpha dra', filename, 51.4167, 2459713.479468, 0.05)

    # Run plot radial velocity for alpha dra
    #filename = "D:\\ASTRO\\Starex\\alphadra\\radial2.lst"
    #plotRadialVelocityCurveSB1(filename, 51.440, -13.5, -47.48, 0.426, -21.80, 0.115)

    #

    # find spec files 
    specs = []
    wdir = 'D:\\ASTRO\\Starex\\alphadra\\'
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile('vr_(.+)_(\d+)_(\d+)(.*).fit')
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : 0 spectrum file found !')
    else:
        logging.info('\U0001F4C1 %d spectra files found !' % (len(specs)))
        data = getBinSysData(specs, 51.440)
        for key, value in data.items():
            data[key]['rv'] = getRv(value['centroid'])
        plotRadialVelocityCurveSB1_2(data, 51.440, -13.5, -47.48, 0.426, -21.80, 0.1)

        print('----')
        for key, value in data.items():
            print(key, value['phase'], value['centroid'])
    
