import logging
import re
import os
import math
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
from binarystarsolve.binarystarsolve import StarSolve


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
    logging.info('---------------------------------------------')
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
        
        t, p = findNextDateByPhase(period, i, jd0, step)
        if(res != None):
            results.append('\U00002705 Phase = %f : next date=> %s (p=%s) (ok : %s)' % ( round(i,2),t.iso, round(p,2), res))
        else:
            results.append('\U0000274C Phase = %f : next date => %s (p=%s)' % (round(i,2), t.iso, round(p,2)))
    
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
            return (t,p)
        t = t + 1*u.d 

def getPhase(jd0, period, jd):
    return (jd-jd0) / period % 1 

def getBinSysData(specs, period):
    obs = {}
    for s in specs:
        f = fits.open(s)
        header = f[0].header
        jd = header['JD-OBS']
        obs[jd] = {'fits':s, 'date':header['DATE-OBS'], 'centroid':float(header['S_CAL']) }

    jd0 = min(obs.keys())
    for jd in obs:
        phase = getPhase(float(jd0), period, float(jd))
        obs[jd]['phase'] = phase
    return obs

def getRadialVelocity(t, t0, K, e, w, v0):
    w = math.radians(w)
    # Mean anomaly
    M = 2 * np.pi *  ((t - t0)%1)    
    # Eccentric anomaly
    E = utilities.eccentric_anomaly_from_mean(e,M, tolerance=0.00001)   
    # True anomaly
    f = utilities.true_anomaly_from_eccentric(e, E) 
    return (K * (e * np.cos(w) + np.cos(w + f)) + v0 ) 

def getRv(deltalambda, lambda0 = 6562.8):
    c = 299792.458
    return  (c  * ((deltalambda-lambda0)/lambda0)) 


def initPlot():
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'

    fig, ax =  plt.subplots(figsize=(9,6))

    #Add Graph title
    #plt.suptitle(r"$\bf{α}$" + " " + r"$\bf{Dra}$" +" - HD123299 - Radial-velocity measurements as a function of Julian Date \nSkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM" ,fontsize=9, fontweight=0, color='black' )
    plt.suptitle(r"$\bf{α}$" + " " + r"$\bf{Dra}$" +" - HD123299 - Phased radial-velocities - 8 observations collected from April to October 2022\nSkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM"  ,fontsize=9, fontweight=0, color='black' )

    #Add X axis label
    ax.set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)

    #Add Y axis label
    ax.set_ylabel('Radial velocity [km $s^{-1}$]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    
    ax.grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

    return (fig, ax)

def plotRadialVelocityCurve(v0, K, e, w, jd0,color="red", label="", lw=0.5, alpha=1):
    model_x = np.arange(-0.1,1.1, 0.005)
    model_y = list(map(lambda x: getRadialVelocity(x,jd0,K,e,w,v0), model_x))
    plt.plot(model_x, model_y, color, alpha=alpha, lw=lw, label=label)

def plotRadialVelocityDotsFromData(specs, color):
    xp = []
    xjd = []
    yrv = []
   
    for jd, s in specs.items():
        xjd.append(jd)
        xp.append(s['phase'])
        yrv.append(s['rv']) 
    
    xp.append(1)
    xjd.append(min(xjd))
    yrv.append(specs[min(xjd)]['rv'])

    plt.plot(xp, yrv, color, lw=0.8)
    #plt.errorbar(xp, yrv,yerr = errors, fmt ='o', color='k', lw=0.2)
 
   
def saveAndShowPlot():
    #plt.legend() 

    plt.tight_layout(pad=1, w_pad=0, h_pad=0)

    plt.xticks(np.arange(-0.1, 1.1, 0.1))
    plt.yticks(np.arange(-50, 60, 10))
    plt.savefig('sandbox/alphadra/hd123299-phased-radial-velocities3.png', dpi=300)
    plt.show()  
      

if __name__ == '__main__':
    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 Planning observations of a binary system - Start \U0001F680')

    # Run for mizar
    #filename = 'D:\\ASTRO\\Starex\\mizar\\time2.lst'
    #filename = '/Volumes/Samsung_T5/ASTRO/Starex/mizar/time2.lst'
    #res = binarySystemObservation('mizar', filename, 20.53835, 2459720.381331, 0.045)

    # Run for alpha dra
    #filename = 'D:\\ASTRO\\Starex\\alphadra\\time2.lst'
    #filename = '/Volumes/Samsung_T5/ASTRO/Starex/alphadra/time2.lst'
    #res = binarySystemObservation('alpha dra', filename, 51.4167, 2459713.479468, 0.05)

    # Run plot radial velocity for alpha dra
    #filename = "D:\\ASTRO\\Starex\\alphadra\\radial2.lst"
    #plotRadialVelocityCurveSB1(filename, 51.440, -13.5, -47.48, 0.426, -21.80, 0.115)

    # Run plot radial velocity for alpha dra
    #filename = 'sandbox/alphadra/radial3.lst'
    #res = binarySystemObservation('beta aur', filename, 3.9600467300, 2453827.19569, 0.1)

    #

    # find spec files 
    specs = []
    wdir = 'D:\\ASTRO\\Starex\\alphadra\\'
    wdir = '/Volumes/Samsung_T5/ASTRO/starex/studies/mizar/'
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile('ok_(.+)_(\d+)_(\d+)(.*).fit')
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : 0 spectrum file found !')
    else:
        logging.info('\U0001F4C1 %d spectra files found !' % (len(specs)))
        P = 51.41891
        data = getBinSysData(specs, P)
        for key, value in data.items():
            data[key]['rv'] = getRv(value['centroid'])
    
        print('---- start plotRadialVelocityCurve --- ')
        for key, value in data.items():
            print('jd : %s  phase : %s  centroid : %s   rv : %s  date:' % (key, round(value['phase'],3), round(value['centroid'],3), round(value['rv'],3)), value['date'])
    
        # me = (9999, 0)
        # for pg in np.arange(51.2, 53.6, 0.01):
        #     params, err, cov = StarSolve(data_file = "sandbox/alphadra/myRVdata.txt", star = "primary", Pguess= pg, covariance = True, graphs=False)
        #     print('p: %s error:%s'%(pg, np.array(err[:4]).sum()))
        #     if(np.array(err[:4]).sum()<me[0]):
        #         me = (np.array(err[:4]).sum(), pg, params, err, cov)

        # print(me)

        #initPlot()
        # (v0, K, e, w, jd0,color="red", label="", lw=0.5, alpha=1)
        
        #[γ, K, ω, e, T0, P, a, f(M)]
        # params, err, cov = StarSolve(data_file = "sandbox/alphadra/myRVdata.txt", star = "primary", Period= P, covariance = True, graphs=False)
        # plotRadialVelocityCurve(-13.5, 47.48, 0.426, 21.80, 0.135, 'r--',  'R. Bischoff, et al.  - Jul, 2017 ' ,1)
        # plotRadialVelocityCurve(params[0], params[1], params[3], params[2], 0.133, 'k-', 'G. Bertrand      - Oct, 2022 ' , 0.8, 0.8)
        # plotRadialVelocityDotsFromData(data, 'ko')
        # saveAndShowPlot()
        # print('[γ, K, ω, e, T0, P, a, f(M)]')
        # print(params)
        # print(err)

        # params, err, cov = StarSolve(data_file = "sandbox/alphadra/myRVdata.txt", star = "primary", Pguess= P, covariance = True, graphs=True)
        # plt.show()
    
    