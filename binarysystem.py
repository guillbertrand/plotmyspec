import logging
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from orbital import utilities

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
        
        if(res != None):
            results.append('\U00002705 Phase = %f : ok => %s' % (round(i,2), res))
        else:
            t = findNextDateByPhase(period, i, jd0, step)
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
        t = t + 1*u.d 
        p =  (t.jd-jd0) % period / period
        if(p > (phase - step) and p < phase + step):
            return t

def getRadialVelocity(t, t0, P, K, e, w, v0):
    # Mean anomaly
    M = 2 * np.pi *  ((t - t0)%1) / 1   
    # Eccentric anomaly
    E = utilities.eccentric_anomaly_from_mean(e,M, tolerance=0.00001)   
    # True anomaly
    f = utilities.true_anomaly_from_eccentric(e, E) 
    return (K * (e * np.cos(w) + np.cos(w + f)) + v0 ) 

def plotRadialVelocityCurveSB1(filename, P, v0, K, e, w):
    xp = []
    xjd = []
    yrv = []
    ''' with open(filename, newline='\n') as f:
        for line in f:
            p = line.rstrip().split(' ')
            xjd.append(float(p[0]))
            xp.append(float(p[1]))
            lambda0 = 6562.82 
            deltalambda = float(p[2])
            rv =  (299792.458 * ((deltalambda/lambda0) ** 2 - 1) / ((deltalambda/lambda0) ** 2 + 1)) 
            yrv.append(rv) '''
    
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'monospace'

    model_x = np.arange(-0.065,1.01, 0.005)
    model_y = list(map(lambda x: getRadialVelocity(x,-0.065,P,K,e,w,v0), model_x))
    fig, ax =  plt.subplots(figsize=(11,6))

     #Add X axis label
    ax.set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)

    #Add Y axis label
    ax.set_ylabel('Radial velocity [km/s]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    
    #plt.plot(xjd, yrv, 'ro')
    plt.plot(model_x, model_y, alpha=1, color="blue", lw=0.6, label="model")

    ax.grid(color='grey', alpha=0.4, linestyle='-', linewidth=0.5, axis='both')
    plt.legend() 

    plt.xticks(np.arange(0, 1.1, 0.1))
    plt.yticks(np.arange(-50, 60, 10))

    plt.show()  

if __name__ == '__main__':
    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 Planning observations of a binary system - Start \U0001F680')

    # Run for mizar
    #filename = '/Volumes/Samsung_T5/ASTRO/Starex/mizar/time2.lst'
    #res = binarySystemObservation('mizar', filename, 20.53835, 2459720.381331, 0.045)

    # Run for alpha dra
    #filename = '/Volumes/Samsung_T5/ASTRO/Starex/alphadra/time2.lst'
    #res = binarySystemObservation('alpha dra', filename, 51.4167, 2459713.479468, 0.1)

    # Run plot radial velocity for alpha dra
    filename = '/Volumes/Samsung_T5/ASTRO/Starex/alphadra/radial.lst'
    plotRadialVelocityCurveSB1(filename, 51.440, -13.5, -47.48, 0.426, -21.80)
    
