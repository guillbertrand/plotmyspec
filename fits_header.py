from astropy.io import fits
import astropy.units as u
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import glob, os
os.chdir("D:\ASTRO\starex\prism")
for file in glob.glob("[!_]*-*.fits"):
    with fits.open(file, mode='update') as hdul:
        hdr = hdul[0].header
        #objname, number = file.split('-')
        objname = hdr['OBJNAME']
        if objname:
            print(file+' >> '+objname)
            # if 'RA' in hdr and 'DEC' in hdr:
            #     hdr['CRVAL1'] = (hdr['RA'], 'approx coord. in RA')
            #     hdr['CRVAL2'] = (hdr['DEC'], 'approx coord. in DEC')
            # else:
            result_table = Simbad.query_object(objname)
            if result_table:
                ra = result_table[0]['RA']
                dec = result_table[0]['DEC']
                c = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))
                crval1, crval2 = c.to_string().split(' ')
                print('Coord ok : ', crval1, crval2)
                hdr['CRVAL1'] = (float(crval1), 'object coord. in RA')
                hdr['CRVAL2'] = (float(crval2), 'object coord. in DEC')
    print('---')


            