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

__version__ = 2.0

''' 
PlotMySpec 
----------
Small script to plot fits spectra using Python and Matplotlib
----------
'''
class PlotMySpec():
    _spectrum_title = ''
    _spectrum_subtitle = ''
    _spectra_path = []
    _spectums_collection = {}
    _compare_mode = False
    _crop = []
    _conf = {}
    _offset = [0]

    def __init__(self, paths, conf):
        self._conf = conf
        self._compare_mode = self._conf["compare_mode"]
        if("crop" in self._conf and self._conf["crop"]):
            self._crop = np.array(self._conf["crop"].split(',')).astype(np.float64)
        if("compare_mode_y_offset" in self._conf and self._conf["compare_mode_y_offset"] != None):
            if(str(self._conf["compare_mode_y_offset"]).__contains__(';')):
                self._offset = self._conf["compare_mode_y_offset"].split(';')
            else:
                self._offset = [self._conf["compare_mode_y_offset"]]

        i = 0
        for path in paths:
            f = fits.open(path)
            head_tail = os.path.split(path)
            spectrum_data = {}
            spectrum_data["filename"] = os.path.splitext(path)[0]
            spectrum_data["basename"] = head_tail[1]
            spectrum_data["header"] = f[0].header
            if not('CRPIX1' in spectrum_data["header"]):
                logging.info('\U0001F5A5 \U0000274C Unable to process %s' % (head_tail[1]))
                continue
            logging.info('\U0001F5A5 \U00002705 Process %s' % (head_tail[1]))
            # Get first pixel reference
            xRef = spectrum_data["header"]['CRPIX1'] - 1
            #Get length of data axis1
            xlength = spectrum_data["header"]['NAXIS1']
            #Get Wavelength pixel step (Angstrom) and convert
            lambdaStep = spectrum_data["header"]['CDELT1']
            self._lambda1 = spectrum_data["header"]['CRVAL1'] - lambdaStep * xRef
            self._lambda2 = spectrum_data["header"]['CRVAL1'] + lambdaStep * (xlength - xRef)
            #Format dataset into astropy quantities
            if(self._offset and len(self._offset)>1):
                flux= (f[0].data + float(self._offset[i]) * (i)) * u.Jy 
            else:
                flux= (f[0].data + float(self._offset[0]) * (i)) * u.Jy 
            wavelength = np.arange(self._lambda1, self._lambda2, lambdaStep) * u.AA 
            # Spectrum construction
            spectrum_data["spec1d"] = Spectrum1D(spectral_axis=wavelength, flux=flux)
            spectrum_data["header"]["DATE-OBS"] = spectrum_data["header"]['DATE-OBS'].split('.')[0]
            
            self._spectums_collection[spectrum_data["header"]["DATE-OBS"]+str(i)] = spectrum_data
            i+=1

    def parsePattern(self, spec, pattern):
        for header_key in spec['header']:
            pattern = pattern.replace('%%'+header_key+'%%', str(spec['header'][header_key]))
        return pattern

    def plotSpec(self):
        for spec in self._spectums_collection.values():
            plt, ax = self.initPlot(spec)
            pngFilename = spec['filename']+'_plot.png'
            ax.plot(spec["spec1d"].spectral_axis, spec["spec1d"].flux, label=spec["header"]['OBJNAME'], alpha=1, color="black", lw=self._conf['line_width']) 
            plt.tight_layout(pad=1, w_pad=0, h_pad=0)
            dpi = self._conf['dpi'] if 'dpi' in self._conf else 150
            plt.savefig(pngFilename, dpi=dpi)
            logging.info('\U0001F4C8 Plot %s > save as %s' % (spec["basename"], pngFilename))
            plt.show()
            if('lines' in self._conf and self._conf['lines']):
                plt, ax = self.initPlot(spec, True)
                pngFilename = spec['filename']+'_plot_wl.png'
                ax.plot(spec["spec1d"].spectral_axis, spec["spec1d"].flux, label=spec["header"]['OBJNAME'], alpha=1, color="black", lw=self._conf['line_width']) 
                plt.tight_layout(pad=1, w_pad=0, h_pad=0)
                dpi = self._conf['dpi'] if 'dpi' in self._conf else 150
                plt.savefig(pngFilename, dpi=dpi)
                logging.info('\U0001F4C8 Plot %s > save as %s' % (spec["basename"], pngFilename))
                plt.show()
            
    def plotSpecGroupMode(self):
        shift = 0
        if('shift' in self._conf):
            shift = self._conf['shift']
        items = list(self._spectums_collection.values())
        plt, ax = self.initPlot(items[0], 'lines' in self._conf)
        pngFilename = items[0]['filename']+'_group_plot.png'
        for key, spec in sorted(self._spectums_collection.items()):
            label = self.parsePattern(spec, self._conf['label_pattern'])
            c = self._conf["compare_mode_color"] if "compare_mode_color" in self._conf and self._conf["compare_mode_color"] else None
            ax.plot(spec["spec1d"].spectral_axis+shift* u.AA, spec["spec1d"].flux, label=label, color=c, alpha=1, lw=self._conf['line_width'])

        if(not "compare_mode_no_label" in self._conf or self._conf["compare_mode_no_label"] != 1):
            if(c and self._offset):
                count = 0
                for key, spec in sorted(self._spectums_collection.items()):
                    label = self.parsePattern(spec, self._conf['label_pattern'])
                    if(self._crop[1]):
                        ax.text(self._crop[1] - 12, 1.07+self._offset[0]*count,label, size='medium')
                    else:
                        ax.text(self._lambda2 - 12, 1.07+self._offset[0]*count,label, size='medium')
                    count += 1
            else:    
                plt.legend() 
       

        


        dpi = self._conf['dpi'] if 'dpi' in self._conf else 150
        plt.savefig(pngFilename, dpi=dpi)
        logging.info('\U0001F4C8 Plot spectra > save as %s' % (pngFilename))
        plt.show()

    def initPlot(self, spec, withlines=False):
        plt.rcParams['font.size'] = self._conf["font_size"]
        plt.rcParams['font.family'] = self._conf["font_family"]

        fig, ax = plt.subplots(figsize=(self._conf["fig_size_x"],self._conf["fig_size_y"]))
       
        self._spectrum_title = ''
        obj = self._conf['object_name'] if(self._conf['object_name']) else spec['header']['OBJNAME']
        split_oname = obj.split(' ')
        for w in split_oname:
            self._spectrum_title += r"$\bf{%s}$ " % (w)
        self._spectrum_title += self.parsePattern(spec, self._conf["title_pattern"])
        self._spectrum_subtitle = self.parsePattern(spec, self._conf["subtitle_pattern"])
        
        #Add Graph title
        plt.title(self._spectrum_title+'\n'+self._spectrum_subtitle,fontsize=self._conf["title_font_size"], fontweight=0, color='black', fontname = self._conf["font_family"])

        #Add X axis label
        ax.set_xlabel(self._conf['x_label'], fontdict=None, labelpad=None, fontname = self._conf["font_family"],size=self._conf["font_size"])

        #Add Y axis label
        ax.set_ylabel(self._conf['y_label'], fontdict=None, labelpad=None, fontname = self._conf["font_family"],size=self._conf["font_size"])
        
        #Add grid 
        if not (self._conf['no_grid']):
            ax.grid(color='grey', alpha=0.4, linestyle='-', linewidth=0.5, axis='both')

        if(withlines and self._conf['lines']):
            for line in self._conf['lines']:
                n = ''
                name, lam, offset_x, offset_y = line.split(',')
                split_oname = name.split(' ')
                for w in split_oname:
                    n += r"$\bf{%s}$ " % (w)
                plt.axvline(x=float(lam), color=self._conf['lines_color'], linestyle='--', linewidth=0.7, alpha=0.6)
                if(self._conf['lines_display_name']):
                    ax.text(float(lam)+float(offset_x), float(offset_y), n+'-'+lam+' Å', rotation=90, va='bottom', color=self._conf['lines_color'],size='7')

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
        logging.info('\U0001F4C1 Error : pms.config.yaml not found !')

    # find spec files 
    wdir = path
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile(conf["spec_file_regex"])
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : 0 spectrum file found !')
    else:
        logging.info('\U0001F4C1 %d spectra files found !' % (len(specs)))
        # Run PlotMySpec
        smp = PlotMySpec(specs, conf)
        smp.run()