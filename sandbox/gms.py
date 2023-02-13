import glob
from PIL import Image

# filepaths
fp_in = "/Volumes/Samsung_T5/ASTRO/Starex/alphadra/ok_hd*_plot_wl.png"
fp_out = "/Volumes/Samsung_T5/ASTRO/Starex/alphadra/anim.gif"

# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
imgs = (Image.open(f) for f in sorted(glob.glob(fp_in)))
img = next(imgs)  # extract first image from iterator
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=450, loop=1)
