from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
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

X,Y,CN,ANALYTICAL,CNERROR = genfromtxt(r'dataCrankN.txt', unpack=True) 
X = X.reshape((mMinus0ne,mMinus0ne))
Y = Y.reshape((mMinus0ne,mMinus0ne))
CN = CN.reshape((mMinus0ne,mMinus0ne))
ANALYTICAL = ANALYTICAL.reshape((mMinus0ne,mMinus0ne))
CNERROR = CNERROR.reshape((mMinus0ne,mMinus0ne))
X=transpose(X)
Y=transpose(Y)
CN=transpose(CN)
ANALYTICAL=transpose(ANALYTICAL)
CNERROR=transpose(CNERROR)

# ------ Plot the surface.
fig = plt.figure()
ax = fig.gca(projection='3d')
surf1 = ax.plot_surface(X, Y, ANALYTICAL, cmap=cm.coolwarm,
                      linewidth=0, antialiased=False,vmin=0,vmax=4.2)
surf2 = ax.plot_surface(X, Y, BACKWARD, cmap=cm.coolwarm,
                      linewidth=0, antialiased=False,vmin=0,vmax=4.2)               
# ------ Customize the z axis.
#ax.set_zlim(300, 400)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$\mathrm{Percentage \ Error}$',labelpad=12,rotation=90)
ax.zaxis.set_rotate_label(False) 
# ------- Add a color bar which maps values to colors.
fig.colorbar(surf1, shrink=0.5, aspect=5)
#fig.colorbar(surf2, shrink=0.5, aspect=5)
# ------- Figure show-save
plt.show()
#plt.tight_layout()
#plt.savefig('/Users/ramytanios/Desktop/Computational Fluid Dynamics Codes/Fortran Codes/Transient Term/results/addSmart.pdf')

