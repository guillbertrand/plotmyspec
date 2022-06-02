# PlotMySpec
Small script to __plot fits spectrums__ using Python and Matplotlib

![11cam spectrum](https://guillaumebertrand.notion.site/image/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2Fb0ab911c-ebf2-485b-90cc-7abda355c68b%2F_11cam_20220516_92_plot.png?table=block&id=8c44c4d1-9b7f-418f-b6b9-56ad589a4f26&spaceId=7d247eda-d75c-46b1-bab6-a26d366d8605&width=2000&userId=&cache=v2)

## Requirements 
```bash
pip install -r requirements.txt
```

## Quickstart

Create a configuration file (pms.yaml) in your working directory with this command
```bash
$ python pms.py 2022-05-14 -i
or 
$ python pms.py 2022-05-14 --init
```

Customize your configuration file (pms.yaml).
__title_pattern, label_pattern and subtitle_pattern__ properties allow you to auto-fill value with your fits file header contents.

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
object_name: '' # by default equal to %%OBJNAME%% prop in your fits fle
title_pattern: '- %%DATE-OBS%% - %%EXPTIME2%% - R=%%SPE_RPOW%% - %OBSERVER%'
label_pattern: '%%DATE-OBS%%'
subtitle_pattern: '%%BSS_INST%%'
spec_file_regex: _(.+)_(\d+)_(\d+)(.*).fit 
crop: 6540,6585,0.25,1.1 # optional
```

## Plot one or more spectrums on distincts graphs

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

## Plot multiple spectrums on the same graph

```bash
$ python pms.py 2022-05-14 -g
or
$ python pms.py 2022-05-14 --group

* INFO - 🚀  PlotMySpec 0.3 - Start 🚀
* INFO - 🖥 Process 2022-05-14/_sheliak_20220513_912 fits file
* INFO - 🖥 Process 2022-05-14/_sheliak_20220513_979 fits file
* INFO - 📈Plot spectrums > save as 2022-05-14/_sheliak_20220513_912_group_hd_plot.png
```

![multiple spectrums](http://www.astrosurf.com/uploads/monthly_2022_05/_sheliak_20220520_956_group_plot.png.2991b5a388ae1a37891d57211ca967dc.png)
