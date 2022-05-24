# PlotMySpec
Small script to __plot fits spectrums__ using Python and Matplotlib


## Requirements 
```bash
pip install -r requirements.txt
```

## Quickstart

Create a configuration file in your working directory with this command
```bash
$ python pms.py 2022-01-14 -i
or 
$ python pms.py 2022-01-14 --init
```

Customize your configuration file.
__*_pattern__ properties allow you to auto-fill value with your fits file header contents.

```bash
---
line_width: 0.7
font_size: 8
title_font_size: 9
font_family: monospace
fig_size_x: 11
fig_size_y: 6
x_label: Wavelength in Å
y_label: Relative intensity
no_grid: 0
display_object_title: 1
title_pattern: '- %%DATE-OBS%% - %%EXPTIME2%% - R=%%SPE_RPOW%% - %OBSERVER%'
label_pattern: '%%DATE-OBS%%'
subtitle_pattern: '%%BSS_INST%%'
spec_file_regex: _(.+)_(\d+)_(\d+)(.*).fit
```

## Plots one or more spectrums on distincts graphs

```bash
$ python pms.py 2022-01-14 

* INFO - 🚀  PlotMySpec 0.3 - Start 🚀
* INFO - 🖥 Process 2022-05-14/_phecda_20220513_912 fits file
* INFO - 🖥 Process 2022-05-14/_hd123299_20220513_979 fits file
* INFO - 📈Plot /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_phecda_20220513_912 fits file > save as 2022-05-14/_phecda_20220513_912_hd_plot.png
Press [enter] to continue.
* INFO - 📈Plot /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_hd123299_20220513_979 fits file > save as 2022-05-14/_hd123299_20220513_979_hd_plot.png
Press [enter] to continue.
```

## Plots multiple spectrums on the same graph

```bash
$ python pms.py 2022-01-14 -g
or
$ python pms.py 2022-01-14 --group

* INFO - 🚀  PlotMySpec 0.3 - Start 🚀
* INFO - 🖥 Process 2022-05-14/_phecda_20220513_912 fits file
* INFO - 🖥 Process 2022-05-14/_hd123299_20220513_979 fits file
* INFO - 📈Plot spectrums > save as 2022-05-14/_phecda_20220513_912_group_hd_plot.png
```

## Plot & Crop your spectrum

```bash
$ python pms.py 2022-01-14 -c 6550,6660,-1,4
or
$ python pms.py 2022-01-14 --crop 6550,6660,-1,4

* INFO - 🚀  PlotMySpec 0.3 - Start 🚀
* INFO - 🖥 Process 2022-05-14/_phecda_20220513_912 fits file
* INFO - 🖥 Process 2022-05-14/_hd123299_20220513_979 fits file
* INFO - 📈Plot /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_phecda_20220513_912 fits file > save as 2022-05-14/_phecda_20220513_912_hd_plot.png
Press [enter] to continue.
* INFO - 📈Plot /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_hd123299_20220513_979 fits file > save as 2022-05-14/_hd123299_20220513_979_hd_plot.png
Press [enter] to continue.
```
