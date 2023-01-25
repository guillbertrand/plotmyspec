# -*- coding: utf-8 -*-
"""
Lucie Leboulleux code
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from skimage.draw import line_aa
from pylab import *
import cv2


path = "sandbox/spikes/"


def circle(img=None, center=(0, 0), rad=10, border=4, color=[1], smooth=2):
    """
        A method to create a circle on a give image.
        img: Expects numpy ndarray of image. 
        center: center of a circle
        rad: radius of a circle
        border: border of the circle, if -ve, circle is filled
        color: color for circle
        smooth: how smooth should our circle be?(smooth * 360 angles in 0 to 360)
    """
    if type(img) == None:
        raise ValueError("Image can not be None. Provide numpy array instead.")
    ix = center[0]+rad
    angles = 360
    cvalue = np.array(color)
    if type(img) != type(None):
        shape = img.shape
        if len(shape) == 3:
            row, col, channels = shape
        else:
            row, col = shape
            channels = 1
        angles = np.linspace(0, 360, 360*smooth)
        for i in angles:
            a = i*np.pi/180
            y = center[1]+rad*np.sin(a) # it is p=h*sin(theta)
            x = center[0]+rad*np.cos(a)
            
            # since we are wroking on image, coordinate starts from (0, 0) onwards and we have to ignore -ve values
            
        
            if border >= 0:
                b = int(np.ceil(border/2))
                
                x1 = np.clip(x-b, 0, shape[0]).astype(np.int32)
                y1 = np.clip(y-b, 0, shape[1]).astype(np.int32)
                x2 = np.clip(x+b, 0, shape[0]).astype(np.int32)
                y2 = np.clip(y+b, 0, shape[1]).astype(np.int32)
            
                
                img[x1:x2, y1:y2] = cvalue

            else:
                x = np.clip(x, 0, shape[0])
                y = np.clip(y, 0, shape[1])
                r, c = int(x), int(y)
                
                if i > 270:
                    img[center[0]:r, c:center[1]] = cvalue
                elif i > 180:
                    img[r:center[0], c:center[1]] = cvalue
                elif i > 90:
                    img[r:center[0], center[1]:c] = cvalue
                elif i > 0:
                    img[center[0]:r, center[1]:c] = cvalue

        return img

def rectangle(img, pt1, pt2, border=2, color=[0]):
    """
        img: Input image where we want to draw rectangle:
        pt1: top left point (y, x)
        pt2: bottom right point
        border: border of line
        color: color of rectangle line,
        returns new image with rectangle.
        
    """
    p1 = pt1
    pt1 = (p1[1], p1[0])
    p2 = pt2
    pt2 = (p2[1], p2[0])
    b = int(np.ceil(border/2))
    cvalue = np.array(color)
    if border >= 0:
        # get x coordinates for each line(top, bottom) of each side
        # if -ve coordinates comes, then make that 0
        x11 = np.clip(pt1[0]-b, 0, pt2[0])
        x12 = np.clip(pt1[0]+b+1, 0, pt2[0])
        x21 = pt2[0]-b
        x22 = pt2[0]+b+1

        y11 = np.clip(pt1[1]-b, 0, pt2[1])            
        y12 = np.clip(pt1[1]+b+1, 0, pt2[1])   
        y21 = pt2[1]-b
        y22 = pt2[1]+b+1
        # right line
        img[x11:x22, y11:y12] = cvalue
        #left line
        img[x11:x22, y21:y22] = cvalue
        # top line
        img[x11:x12, y11:y22] = cvalue
        # bottom line
        img[x21:x22, y11:y22] = cvalue
        
    else:
        pt1 = np.clip(pt1, 0, pt2)
        img[pt1[0]:pt2[0]+1, pt1[1]:pt2[1]+1] = cvalue
        
    return img

sz = 400
rad = 100
base = np.zeros((sz,sz))

figures = []

c = circle(base.copy(), (200,200), 50, -2, [1])
figures.append(["Refractor", c])

n0 = circle(base.copy(), (200,200), 100, -2, [1])
n0 = rectangle(n0, (198,0), (202,200),-1, [0])


rr, cc, val = line_aa(199, 200, 398, 399)
n0[rr, cc] = val * 0
rr, cc, val = line_aa(200, 199, 399, 398)
n0[rr, cc] = val * 0

rr, cc, val = line_aa(200, 200, 398, 0)
n0[rr, cc] = val * 0
rr, cc, val = line_aa(201, 201, 399, 1)
n0[rr, cc] = val * 0
n1 = circle(n0, (200,201), 17, -2, [0])
figures.append(["Newton 114 f/8", n0])

n1 = circle(base.copy(), (200,200), 100, -2, [1])
n1 = rectangle(n1, (197,0), (201,198),-1, [0])
n1 = rectangle(n1, (0,194), (400,198),-1, [0])
n1 = circle(n1, (198,200), 22, -2, [0])
figures.append(["Newton 400 f/4.3", n1])

n2 = circle(base.copy(), (200,200), 100, -2, [1])
n2 = rectangle(n2, (198,0), (202,400),-1, [0])
n2 = rectangle(n2, (0,198), (400,202),-1, [0])
n1 = circle(n2, (200,201),25, -2, [0])
figures.append(["Newton 500 f/4", n2])

def getPSF(figure):
    Eref = np.fft.fft2(figure) 
    Eref = np.fft.fftshift(Eref)
    psf = np.abs(Eref)**2
    norm = np.max(psf)
    psf=psf/norm
    psfLog = np.log10(psf)
    return psfLog

plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'monospace'
fig, ax = plt.subplots(len(figures), 2, figsize=(4, 8))
for idx, (title, f) in enumerate(figures):
    psf = getPSF(f)
    ax[idx][0].imshow(f, cmap = cm.Greys_r)
    ax[idx][0].axis('off')
    ax[idx][1].imshow(psf, cmap = cm.inferno,vmin=-4.5,vmax=0.)
    ax[idx][1].axis('off')
    ax[idx][0].set_title(title)

plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=1)
plt.savefig(path+'spikes.png', dpi=300, bbox_inches='tight')





