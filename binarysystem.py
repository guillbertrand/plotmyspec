import logging
from astropy.time import Time
import astropy.units as u
from numpy import arange
from datetime import date


def binarySystemObservation(obj, filename, period, jd0, step=0.05):
    logging.info('\U00002728 Start scan for %s period (%s days)' % (obj, period))
    data = []
    results = []
    with open(filename, newline='\n') as f:
        for line in f:
            p = float(line.rstrip().split(' ')[1])
            if (p>=0):
                data.append(p)
        
    for i in arange(0,1,step):
        res = None
        for d in data:
            if d > i-step and d < i+step:
                res = d

        if(res):
            results.append('\U00002705 Phase = %s : ok => %s' % (round(i,2), res))
        else:
            t = findNextDateByPhase(period, i, jd0, step)
            results.append('\U0000274C Phase = %s : next date => %s' % (round(i,2), t.iso))
    
    for r in results:
        logging.info(r)
            
def findNextDateByPhase(period, phase, jd0, step):
    today = date.today()
    time = '%s-%s-%sT%s:%s:00' % (today.strftime("%Y"), today.strftime("%m"), today.strftime("%d"), 23, 00)
    t = Time(time, format='isot', scale='utc')
    res = None
    while not res:
        t = t + 1*u.d 
        p =  (t.jd-jd0) % period / period
        if(p > phase - step and p < phase + step):
            return t

if __name__ == '__main__':

    # Path of your time.lst file from ISIS
    # example time.lst : 
    # @m9 -0.318175
    # @m1 0
    # @m2 0.194897
    # @m3 0.340229
    # @m4 0.392226
    # @m5 0.438418
    # @m6 0.485914
    # @m7 0.534944
    # @m8 0.633233
    # @m9 0.681825
    # @m1 1
    

    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 Planning observations of a binary system - Start \U0001F680')

    # Run for mizar
    filename = '/Volumes/Samsung_T5/ASTRO/Starex/mizar/time2.lst'
    res = binarySystemObservation('mizar', filename, 20.53835, 2459720.381331, 0.03)

    # Run for alpha dra
    #filename = '/Volumes/Samsung_T5/ASTRO/Starex/alphadra/time2.lst'
    #res = completeBinSystemObservation('alpha dra', filename, 51.4167, 2459713.479468, 3)

 
