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
mMinus0ne = int(genfromtxt(r'mMinus0ne.txt'))
X,Y,Z = genfromtxt(r'test.txt', unpack=True,delimiter=',')
#[X,Y] = meshgrid(X,Y) 
X = X.reshape((mMinus0ne,mMinus0ne))
Z = Z.reshape((mMinus0ne,mMinus0ne))
Y = Y.reshape((mMinus0ne,mMinus0ne))
Z=transpose(Z)
X=transpose(X)
Y=transpose(Y)


# ------ Plotting
f = plt.pcolormesh(X,Y,Z,cmap='jet',shading='gouraud')#origin='lower',cmap='jet',interpolation='none')
#f.set_clim([0, 1])
plt.colorbar(f, orientation = 'horizontal',fraction=0.045, pad=0.1)


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
#plt.savefig('/Users/ramytanios/Desktop/Computational Fluid Dynamics Codes/Fortran Codes/Diffusion Convection/results/addSmart.pdf')
