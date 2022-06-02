import os
import re
import yaml
import logging
import argparse

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from specutils import Spectrum1D

__version__ = 1.1


''' 
PlotMySpec 
----------
Small script to plot fits spectrums using Python and Matplotlib
----------
'''
class PlotMySpec():
    _spectrum_title = ''
    _spectrum_subtitle = ''
    _spectrums_path = []
    _spectums_collection = []
    _compare_mode = False
    _crop = []
    _conf = {}

    def __init__(self, paths, conf):
        self._conf = conf
        self._compare_mode = self._conf["compare_mode"]
        if("crop" in self._conf and self._conf["crop"]):
            self._crop = np.array(self._conf["crop"].split(',')).astype(np.float64)

        for path in paths:
            f = fits.open(path)
            head_tail = os.path.split(path)
            spectrum_data = {}
            spectrum_data["filename"] = os.path.splitext(path)[0]
            spectrum_data["header"] = f[0].header
            if not('CRPIX1' in spectrum_data["header"]):
                continue
            logging.info('\U0001F5A5 Process %s fits file' % (head_tail[1]))
            # Get first pixel reference
            xRef = spectrum_data["header"]['CRPIX1'] - 1
            #Get length of data axis1
            xlength = spectrum_data["header"]['NAXIS1']
            #Get Wavelength pixel step (Angstrom) and convert
            lambdaStep = spectrum_data["header"]['CDELT1']
            self._lambda1 = spectrum_data["header"]['CRVAL1'] - lambdaStep * xRef
            self._lambda2 = spectrum_data["header"]['CRVAL1'] + lambdaStep * (xlength - xRef)
            #Format dataset into astropy quantities
            flux= f[0].data * u.Jy
            wavelength = np.arange(self._lambda1, self._lambda2, lambdaStep) * u.AA
            # Spectrum construction
            spectrum_data["spec1d"] = Spectrum1D(spectral_axis=wavelength, flux=flux)
            spectrum_data["header"]["DATE-OBS"] = spectrum_data["header"]['DATE-OBS'].split('.')[0]
            
            self._spectums_collection.append(spectrum_data)

    def parsePattern(self, spec, pattern):
        for header_key in spec['header']:
            pattern = pattern.replace('%%'+header_key+'%%', str(spec['header'][header_key]))
        return pattern

    def plotSpec(self):
        for spec in self._spectums_collection:
            plt, ax = self.initPlot(spec)
            pngFilename = spec['filename']+'_hd_plot.png'
            pngLRFilename = spec['filename']+'_plot.png'
            ax.plot(spec["spec1d"].spectral_axis, spec["spec1d"].flux, label=spec["header"]['OBJNAME'], alpha=1, color="black", lw=self._conf['line_width']) 
            plt.tight_layout(pad=1, w_pad=0, h_pad=0)
            plt.savefig(pngFilename, dpi=300)
            plt.savefig(pngLRFilename, dpi=150)
            logging.info('\U0001F4C8 Plot %s fits file > save as %s' % (spec["filename"], pngFilename))
            plt.show()
            
    def plotSpecGroupMode(self):
        plt, ax = self.initPlot(self._spectums_collection[0])
        pngFilename = self._spectums_collection[0]['filename']+'_group_hd_plot.png'
        pngLRFilename = self._spectums_collection[0]['filename']+'_group_plot.png'
        for spec in self._spectums_collection:
            label = self.parsePattern(spec, self._conf['label_pattern'])
            ax.plot(spec["spec1d"].spectral_axis, spec["spec1d"].flux, label=label, alpha=1, lw=self._conf['line_width']) 
        plt.legend() 
        plt.savefig(pngFilename, dpi=300)
        plt.savefig(pngLRFilename, dpi=150)
        logging.info('\U0001F4C8 Plot spectrums > save as %s' % (pngFilename))
        plt.show()

    def initPlot(self, spec):
        plt.rcParams['font.size'] = self._conf["font_size"]
        plt.rcParams['font.family'] = self._conf["font_family"]

        fig, ax = plt.subplots(figsize=(self._conf["fig_size_x"],self._conf["fig_size_y"]))
       
        obj = self._conf['object_name'] if(self._conf['object_name']) else spec['header']['OBJNAME']
        self._spectrum_title = r"$\bf{%s}$ " % (obj)
        self._spectrum_title += self.parsePattern(spec, self._conf["title_pattern"])
        self._spectrum_subtitle = self.parsePattern(spec, self._conf["subtitle_pattern"])
        
        #Add Graph title
        plt.suptitle(self._spectrum_title,fontsize=self._conf["title_font_size"], fontweight=0, color='black', x=0.515,y=0.97, fontname =self._conf["font_family"])
        plt.title(self._spectrum_subtitle,fontsize=self._conf["font_size"], fontweight=0, color='black',y=1.03, fontname = self._conf["font_family"])

        #Add X axis label
        ax.set_xlabel(self._conf['x_label'], fontdict=None, labelpad=None, fontname = self._conf["font_family"],size=self._conf["font_size"])

        #Add Y axis label
        ax.set_ylabel(self._conf['y_label'], fontdict=None, labelpad=None, fontname = self._conf["font_family"],size=self._conf["font_size"])
        
        #Add grid 
        if not (self._conf['no_grid']):
            ax.grid(color='grey', alpha=0.4, linestyle='-', linewidth=0.5, axis='both')
        
        # Crop spectrum 
        ax.set_xlim(self._crop[:2])  if(len(self._crop) >= 2) else ax.set_xlim([self._lambda1, self._lambda2])
        if(len(self._crop) == 4):
            ax.set_ylim(self._crop[2:])
            
        plt.tight_layout(pad=1, w_pad=0, h_pad=0)

        return (plt, ax)

    def run(self):
        if(self._compare_mode and len(self._spectums_collection) > 1):
            self.plotSpecGroupMode()
        else:
            self.plotSpec()
        

if __name__ == '__main__':
    specs = []
    #
    default_conf_filename = 'pms.config.yaml'

    #
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to your fits directory")
    parser.add_argument("-c", "--config", type=str, default=default_conf_filename, help="Custom config file name")
    
    args = parser.parse_args()
    path = args.path
    conf_filename = args.config

    wdir = path

    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 PlotMySpec %s - Start \U0001F680' % __version__)

    # load yaml configuration file
    confpath = os.path.join(wdir,conf_filename)
    cust_confpath = os.path.join(wdir, conf_filename)
    if not (os.path.exists(confpath)):
        confpath = conf_filename
        if not (os.path.exists(confpath)):
            confpath = None

    if confpath:
        logging.info('\U00002728 Load configuration file \U0001F527  %s' % (confpath))
        with open(confpath, 'r', encoding='utf8') as f:
            conf = yaml.load(f,Loader=yaml.FullLoader)
    else :
        logging.info('\U0001F4C1 Error : pms.config.yaml not found !' % (len(specs)))

    # find spec files 
    wdir = path
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile(conf["spec_file_regex"])
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : %d spectrum file found !' % (len(specs)))
    else:
        logging.info('\U0001F4C1 %d spectrums files found !' % (len(specs)))
        # Run PlotMySpec
        smp = PlotMySpec(specs, conf)
        smp.run()