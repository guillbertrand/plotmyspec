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
    return (jd-jd0) % period / period

def getBinSysData(specs, period, jd0):
    obs = {}
    for s in specs:
        f = fits.open(s)
        header = f[0].header
        jd = header['JD-OBS']
        obs[jd] = {'fits':s, 'centroid1':float(header['S_CAL'].split(';')[0]), 'centroid2':float(header['S_CAL'].split(';')[1]) }

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
    return  (c  * (deltalambda/lambda0)) - c


def initPlot():
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'

    fig, ax =  plt.subplots(figsize=(9,6))

    #Add Graph title
    plt.suptitle(r"$\bf{Mizar}$" +" - HD116656 - Phased radial-velocities - %s observations collected from May to June 2022" % len(specs),fontsize=9, fontweight=0, color='black' )
    plt.title("SkyWatcher refractor D=72mm f/6 + Star'Ex (2400 l/mm, 80x125, 10 μm slit) + ASI 183MM",fontsize=8, fontweight=0, color='black')

    #Add X axis label
    ax.set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)

    #Add Y axis label
    ax.set_ylabel('Radial velocity [km/s]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    
    ax.grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both')

    return (fig, ax)

def plotRadialVelocityCurve(v0, K, e, w, jd0,style="-",color="red", label="", lw=0.5, alpha=1):
    model_x = np.arange(0,1.1, 0.005)
    model_y = list(map(lambda x: getRadialVelocity(x,jd0,K,e,w,v0), model_x))
    plt.plot(model_x, model_y, style, alpha=alpha, lw=lw, label=label, color=color)

def plotRadialVelocityDotsFromData(specs, color, period):
    xp = []
    xjd = []
    yrv1 = []
    yrv2 = []
   
    for jd, s in specs.items():
        xjd.append(jd)
        xp.append(s['phase'])
        yrv1.append(s['rv1']) 
        yrv2.append(s['rv2']) 
    
    plt.plot(xp, yrv1, color, lw=0.8)
    plt.plot(xp, yrv2, color, lw=0.8)
 
   
def saveAndShowPlot():
    plt.legend() 

    plt.tight_layout(pad=1, w_pad=0, h_pad=0)

    plt.xticks(np.arange(0, 1.1, 0.1))
    plt.yticks(np.arange(-90, 90, 10))
    plt.savefig('sandbox/mizar/mizar-phased-radial-velocities2.png', dpi=300)
    plt.show()  
      

if __name__ == '__main__':
    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    #

    period = 20.53835
    jd0 = 2459720.381331

    #find spec files 
    specs = []
    wdir = 'D:\\ASTRO\\Starex\\mizar\\'
    wdir = 'sandbox/mizar/'
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile('rv(\d+).fit')
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : 0 spectrum file found !')
    else:
        logging.info('\U0001F4C1 %d spectra files found !' % (len(specs)))
        data = getBinSysData(specs, period, jd0)
        for key, value in data.items():
            data[key]['rv1'] = getRv(value['centroid1'])
            data[key]['rv2'] = getRv(value['centroid2'])

        print('---- start plotRadialVelocityCurve --- ')
        for key, value in data.items():
            print('jd : %s  phase : %s  centroid 1 : %s, centroid 2 : %s, rv1 : %s, rv2 : %s' % (key, round(value['phase'],3), round(value['centroid1'],3), round(value['centroid2'],3), round(value['rv1'],3), round(value['rv2'],3)))

        print('---- start output 1 --- ')
        with open('sandbox/mizar/myRVdata2.txt', 'w') as f: 
            for key, value in data.items():
                if(value['phase']>0.6):
                    output = '%s %s %s' % (key-2400000, round(value['rv2'],3),round(value['rv1'],3))
                else:
                    output = '%s %s %s' % (key-2400000, round(value['rv1'],3),round(value['rv2'],3))
                f.write(output+'\n')
                print(value['fits'], output)


    initPlot()
    v0 = 0.625
    plotRadialVelocityCurve(-5.6, 68.85, 0.542, 104.16, v0, '-', 'red',  'Primary - Budovicová et al. - 2004',1)
    plotRadialVelocityCurve(-5.6, -65.51, 0.542, 104.16, v0, '--', 'red',  'Secondary - Budovicová et al. - 2004',1)
    
    

    from binarystarsolve.binarystarsolve import StarSolve
    params, err = StarSolve(data_file = "sandbox/mizar/myRVdata2.txt", star = "both", Pguess=period,  covariance = False, graphs=False)
  
    # [γ, K, ω, e, T0, P, a, f(M)]
    print(params[0])
    print(params[1])

    plotRadialVelocityCurve(params[1][0], params[1][1], params[1][3], params[1][2], v0+0.005, '--', 'black', 'Primary - G. Bertrand - Jun, 2022', 0.8, 0.8)
    plotRadialVelocityCurve(params[0][0], params[0][1], params[0][3], params[0][2], v0+0.005, '-', 'black',  'Secondary - G. Bertrand - Jun, 2022', 0.8, 0.8)

    
    plotRadialVelocityDotsFromData(data, 'ko', period)
    saveAndShowPlot()

    print(err[0])
    print(err[1])
