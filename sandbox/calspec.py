import numpy as np

from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler

from astropy import units as u 
import astropy.wcs as fitswcs

from getCalspec.getCalspec import *

def resample_calspec(objname, wvl = (3000, 11000), delta= 0.01, filename="spec.fits"):
    c = Calspec(objname)
    n = c.get_spectrum_numpy()  
    s = Spectrum1D(spectral_axis=n["WAVELENGTH"], flux=n["FLUX"])
    resample_grid = np.arange(wvl[0], wvl[1], delta)
    fluxc_resample = LinearInterpolatedResampler()
    output_spectrum1D = fluxc_resample(s, resample_grid* u.AA) 

    my_wcs = fitswcs.WCS(header={'CDELT1': delta, 'CRVAL1': wvl[0], 'CUNIT1': 'Angstrom', 'CRPIX1': 1.})
    resampleSpec = Spectrum1D(flux=output_spectrum1D.flux, wcs=my_wcs)
    resampleSpec.write(filename, overwrite=True, format="wcs1d-fits")

resample_calspec('hd120315', filename="betauma.fits")