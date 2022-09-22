import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy.special import sph_harm
from itertools import product, combinations
from scipy.integrate import simps

import glob
from PIL import Image

#-------------------------------------------------------------------------------
# Longitudinal and surface magnetic field calculation --------------------------
#-------------------------------------------------------------------------------
def Beff(phi,inc,beta,Bd,x0=0,y0=0,z0=0,phase_offset=0,u=0.3):

	#Length of arrays
	PH = len(phi) #Rotational phase
	N = 100 #Spatial grid size

	#Defining polar and azimuthal grids
	MU = np.linspace(0.,1.,N) 
	PHI = np.linspace(0.,1.,N)*2.*np.pi

	#Creating 2D meshgrids
	MU_grid, PHI_grid = np.meshgrid( MU, PHI, indexing='xy')
	X = np.sqrt(1.-MU_grid**2)*np.cos(PHI_grid) 
	Y = np.sqrt(1.-MU_grid**2)*np.sin(PHI_grid) 
	Z = MU_grid 

	#Variable setup
	Bl = np.zeros(PH) #empty array for longitudinal field
	Bs = np.zeros(PH) #empty array for surface field
	Bs2 = []
	inc = np.radians(inc)
	beta = np.radians(beta)
	phi = 2.*np.pi * (phi + phase_offset)
	for ph in range(0,PH):
		
		#Observer's frame (X,Y,Z)
		#Inclined and rotated frame (XR,YR,ZR)
		#Magnetic frame (XM,YM,ZM)

		#Set of transformations: observer's frame -> magnetic frame
		XR=X*np.cos(inc)-Z*np.sin(inc)
		YR=Y
		ZR=X*np.sin(inc)+Z*np.cos(inc)
		XM=(XR*np.cos(phi[ph])+YR*np.sin(phi[ph]))*np.cos(beta)+ZR*np.sin(beta)
		YM=-XR*np.sin(phi[ph])+YR*np.cos(phi[ph])
		ZM=-(XR*np.cos(phi[ph])+YR*np.sin(phi[ph]))*np.sin(beta)+ZR*np.cos(beta)

		#Dipole offset (x0,y0,z0)
		x = XM-x0
		y = YM-y0
		z = ZM-z0
		r = np.sqrt(x**2 + y**2 + z**2)

		#Components of the dipolar magnetic field vector in cartesian coordinates
		Bx = 3.*Bd/2*x*z/r**5
		By = 3.*Bd/2*y*z/r**5
		Bz = Bd/2*( 3.*z**2 - r**2 )/r**5

		#Set of transformations: magnetic frame -> observer's frame
		BXR = (Bx*np.cos(beta) - Bz*np.sin(beta))*np.cos(phi[ph]) - By*np.sin(phi[ph])
		BYR = (Bx*np.cos(beta) - Bz*np.sin(beta))*np.sin(phi[ph]) + By*np.cos(phi[ph])
		BZR = Bx*np.sin(beta) + Bz*np.cos(beta)
		BX = BXR*np.cos(inc) + BZR*np.sin(inc) 
		BY = BYR
		BZ = -BXR*np.sin(inc)  + BZR*np.cos(inc) 

		#Magnitude of the field
		B = np.sqrt(BX**2 + BY**2 + BZ**2)	

		#Surface integrals over the visible portion of the star
		meanBz = simps(simps(BZ*MU_grid*(1-u+u*MU_grid),MU),PHI)
		meanBs = simps(simps(B*MU_grid*(1-u+u*MU_grid),MU),PHI)
		mean = simps(simps(MU_grid*(1-u+u*MU_grid),MU),PHI) 

		#Resulting longitudinal and surface magentic field 
		Bl[ph] = meanBz/mean 
		Bs[ph] = meanBs/mean 

	return [Bl,Bs]


Bd = 3500

q = 250

phi = np.linspace(0, np.pi, q)
theta = np.linspace(0, 2*np.pi, q)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta) 
y = np.sin(phi) * np.sin(theta) 
z = np.cos(phi)

m, l = 1, 1



# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = sph_harm(m, l, theta, phi).real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'monospace'
# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')

def rotateY(x,y,z,a=2.72): #155.84°
    t = np.transpose(np.array([x,y,z]), (1,2,0))
    m = [[math.cos(a), 0, math.sin(a)],[0,1,0],[-math.sin(a), 0, math.cos(a)]]
    x,y,z = np.transpose(np.dot(t, m), (2,0,1))
    return x,y,z


pole = rotateY([[0.]],[[0.]],[[1.]])
# rotate the samples by pi / 4 radians around y
x,y,z = rotateY(x,y,z)

tmpsurf = ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.RdBu(fcolors), shade= False)
tmpsurf.set_sort_zpos(0.0)

tmpLineax1 = ax.quiver(0.0, 0.0, 1.0, 0.0, 0.0, .5, normalize=False, linewidth=2, color=(0,0,0), pivot='tail',arrow_length_ratio=0.0)
tmpLineax2 = ax.quiver(0.0, 0.0,-1.0, .0, 0.0,-0.5, normalize=False, linewidth=2, color=(0,0,0), pivot='tail',arrow_length_ratio=0.0)
#And try to build an equator out of quiver objects too
eqR = np.max(1.01)
for i in range(q):
    theta = 2*np.pi*(i/q)
    dtheta = 2*np.pi*(1./q)
    circ_x1 = eqR*np.sin(theta)
    circ_y1 = eqR*np.cos(theta)
    circ_dx = eqR*(np.sin(theta+dtheta) - np.sin(theta))
    circ_dy = eqR*(np.cos(theta+dtheta) - np.cos(theta))
    tmpCirc = ax.quiver(circ_x1, circ_y1, 0.0, circ_dx, circ_dy, 0.0, normalize=False, linewidth=1, color=(0.0,0.0,0.0), pivot='tail',arrow_length_ratio=0.0)

phase  = 0

# Make axes limits 
xyzlim = np.array([ax.get_xlim3d(),ax.get_ylim3d(),ax.get_zlim3d()]).T
XYZlim = [min(xyzlim[0]),max(xyzlim[1])]
ax.set_xlim3d(XYZlim)
ax.set_ylim3d(XYZlim)
ax.set_zlim3d(XYZlim)

ax.set_box_aspect((1, 1, 1))

#Creates array for colorbar from 0 to 1. 
a = np.array( [-3.5, 3.5])

#Creates colorbar
m = cm.ScalarMappable(cmap=cm.RdBu_r)
m.set_array(a)
plt.colorbar(m)

text = ax.text2D(0.45, 0.05, 'φ {:0.2f}'.format(0), transform=ax.transAxes, fontsize='large')

# Turn off the axis planes
ax.set_axis_off()
plt.tight_layout(pad=1, w_pad=0, h_pad=0)
#plt.suptitle("Distribution of the magnetic field on the surface\n of " +r"$\bf{α2}$" + " " + r"$\bf{CVn}$"+" derived from the mean longitudinal field.\n July 2022, G. Bertrand",fontsize=9, fontweight=0, color='black' )
ax.scatter(pole[0][0]+0.01, pole[1][0]+0.01, pole[2][0]+0.01, color="k", marker="x")


for idx, ph in enumerate(np.linspace(0, 0.99, 100)):
    ax.view_init(90-120, -ph*360)
    #And label this rotation phase
    text.set_text('φ {:0.2f}'.format(ph))
    plt.draw()
    plt.savefig('sandbox/orm/orm4-%04d.png' % idx, dpi=300) 
    break

# filepaths
fp_in = "sandbox/orm/orm4-*.png"
fp_out = "sandbox/orm/anim.gif"

imgs = []
for f in sorted(glob.glob(fp_in)):
    o = Image.open(f)
    o = o.convert('RGB')
    imgs.append(o)
img = imgs[0]  # extract first image from iterator
img.save(fp=fp_out, format='GIF',  append_images=imgs, save_all=True, duration=100, loop=1)
