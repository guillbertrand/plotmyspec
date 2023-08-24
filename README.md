# PlotMySpec
Small script to __plot fits spectra__ using Python and Matplotlib

![11cam spectrum](https://guillaumebertrand.notion.site/image/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2Fb0ab911c-ebf2-485b-90cc-7abda355c68b%2F_11cam_20220516_92_plot.png?table=block&id=8c44c4d1-9b7f-418f-b6b9-56ad589a4f26&spaceId=7d247eda-d75c-46b1-bab6-a26d366d8605&width=2000&userId=&cache=v2)

## Requirements 
```bash
pip install -r requirements.txt
```

## Quickstart

Customize your configuration file (pms.yaml).
__title_pattern, label_pattern and subtitle_pattern__ properties allow you to auto-fill value with your fits file header contents.

```bash
---
# ----------------------------------
# matplotlib general config
# ----------------------------------
line_width: 0.6
font_size: 9
title_font_size: 9
font_family: arial
math_font_family: custom
fig_size_x: 9
fig_size_y: 5
x_label: Wavelength in Å
y_label: Relative intensity
no_grid: 0
dpi: 150

# ----------------------------------
# spectrum config
# ----------------------------------
# object_name:  :auto-filled with keyword OBJNAME from header fits if empty
object_name: 
title_pattern: "- %%DATE-OBS%% - %%EXPTIME2%% - %%OBSERVER%%"
label_pattern: "%%DATE-OBS%%"
subtitle_pattern: "%%BSS_INST%%"
spec_file_regex: '_(.+)_(\d+)_(\d+)(.*).fit'
# crop: allow to crop spectra (ex: "crop: 6500,6600,0,2" --> lambdamin,lambdamax,fluxmin,fluxmax )
crop: 0
# compare_mode : display all spectra on the same plot if active
compare_mode : 0
compare_mode_y_offset : 0
compare_mode_color_cycle: #red,black,yellow... or nothing for automatic mode
compare_mode_no_label: 0

# ----------------------------------
# optional : display vertical lines
# ----------------------------------
lines_display_name: 1
lines_color: red
lines: 
  # line name, position in angstroms, offset x, offset y
  - Hα,6562.82,6.,0.3
```

## Plot one or more spectra on distincts graphs

```bash
$ python pms.py /Volumes/Samsung_T5/2022-01-14/ 

- 🚀 PlotMySpec 1.1 - Start 🚀
- ✨ Load configuration file 🔧  pms.config.yaml
- 📁 4 spectra files found !
- 🖥 ❌ Unable to process _phecda_20220513_912_2D.fits
- 🖥 ✅ Process _phecda_20220513_912.fits
- 🖥 ❌ Unable to process _hd123299_20220513_979_2D.fits
- 🖥 ✅ Process _hd123299_20220513_979.fits
- 📈 Plot _phecda_20220513_912.fits > save as /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_phecda_20220513_912_hd_plot.png
- 📈 Plot _hd123299_20220513_979.fits > save as /Volumes/Samsung_T5/ASTRO/Starex/2022-05-14/_hd123299_20220513_979_hd_plot.png
```

## Plot multiple spectra on the same graph

Set __compare_mode=1__ in your configuration file.

```bash
$ python pms.py /Volumes/Samsung_T5/2022-01-14/ 

- 🚀 PlotMySpec 1.1 - Start 🚀
- ✨ Load configuration file 🔧  /Volumes/Samsung_T5/ASTRO/Starex/mizar/pms.config.yaml
- 📁 7 spectra files found !
- 🖥 ✅ Process @m1.fits
- 🖥 ✅ Process @m2.fits
- 🖥 ✅ Process @m3.fits
- 🖥 ✅ Process @m4.fits
- 🖥 ✅ Process @m5.fits
- 🖥 ✅ Process @m6.fits
- 🖥 ✅ Process @m7.fits
- 📈 Plot spectra > save as /Volumes/Samsung_T5/ASTRO/Starex/mizar/@m1_group_hd_plot.png
```

![multiple spectra](http://www.astrosurf.com/uploads/monthly_2022_05/_sheliak_20220520_956_group_plot.png.2991b5a388ae1a37891d57211ca967dc.png)


![multiple spectra](https://www.notion.so/image/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2F39b32ff5-a445-4bb0-bd0b-a184f4da0c2d%2Fm1_group_plot.png?id=3c2eeaa8-bc58-47e3-9bb0-b5d7d6493cd8&table=block&spaceId=7d247eda-d75c-46b1-bab6-a26d366d8605&width=2000&userId=5d8ddb97-4aa9-4c82-b8c6-1928cd30ab12&cache=v2)

![multiple spectra](https://www.notion.so/image/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2F61475ef3-1ce4-441e-bc55-a0908d90512d%2Fok_hd123299_20220513_979_group_plot.png?id=10488798-79db-4d9f-9af5-6c4b25ba09ef&table=block&spaceId=7d247eda-d75c-46b1-bab6-a26d366d8605&width=2000&userId=5d8ddb97-4aa9-4c82-b8c6-1928cd30ab12&cache=v2)


