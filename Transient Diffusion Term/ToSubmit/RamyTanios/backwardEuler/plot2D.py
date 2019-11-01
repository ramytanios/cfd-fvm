from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc


# ----- Font
rc('font',size=11)
rc('font',family='serif')
rc('axes',labelsize=11)
plt.rcParams["mathtext.fontset"] = "stixsans"


# ----- Make data.
mMinus0ne = 25
X,Y,BACKWARD,ANALYTICAL,BACKERROR = genfromtxt(r'dataBackward.txt', unpack=True) 
X = X.reshape((mMinus0ne,mMinus0ne))
Y = Y.reshape((mMinus0ne,mMinus0ne))
BACKWARD = BACKWARD.reshape((mMinus0ne,mMinus0ne))
ANALYTICAL = ANALYTICAL.reshape((mMinus0ne,mMinus0ne))
BACKERROR = BACKERROR.reshape((mMinus0ne,mMinus0ne))
X=transpose(X)
Y=transpose(Y)
BACKWARD=transpose(BACKWARD)
ANALYTICAL=transpose(ANALYTICAL)
BACKERROR=transpose(BACKERROR)


# ------ Plotting
f = plt.imshow(BACKWARD,origin='lower',cmap='jet',interpolation='bilinear',extent=(0,1,0,1))
#f = plt.contour(Z,origin='lower',cmap='jet',extent=(0,1,0,1))
#plt.clabel(f, inline=1, fontsize=10)
#f.set_clim([0, 1])
plt.colorbar(f, orientation = 'vertical',fraction=0.045, pad=0.1)#format='%.0e')

# ------ Plot stuff on the graph
#plt.quiver([70], [50], [20], [-20], angles='xy', scale_units='xy', scale=1)
#plt.quiver([30], [50], [20], [-20], angles='xy', scale_units='xy', scale=1)


# ------xlabel,ylabel,title,texts
plt.axis('on')
#plt.xlabel(r'$L_x = 1m$')
#plt.ylabel(r'$L_y = 1m$')
#plt.text(30, 102, r'$L_x = 1\mathrm{m}, \ \phi_{in}=1$', fontsize=20)
#plt.text(-8, 65, r'$L_y = 1\mathrm{m}, \ \phi_{in}=0$', fontsize=20,rotation=90)
#plt.text(70, 50, r'$\mathbf{v}$', fontsize=20,rotation=-45)
#plt.title(r'$L_x = 1\mathrm{m}, \ \phi_{in}=1$')


# ------Figure show-save
plt.show()
#plt.tight_layout()
#plt.savefig('/Users/ramytanios/Desktop/Computational Fluid Dynamics Codes/Fortran Codes/Transient Diffusion Term/results/0.1.pdf')
