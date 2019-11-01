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
mMinus0ne = int(genfromtxt(r'mMinus0ne.txt'))
X,Y,Z = genfromtxt(r'test.txt', unpack=True,delimiter=',') 
X = X.reshape((mMinus0ne,mMinus0ne))
Z = Z.reshape((mMinus0ne,mMinus0ne))
Y = Y.reshape((mMinus0ne,mMinus0ne))
Z=transpose(Z)
X=transpose(X)
Y=transpose(Y)

# ------ Plot the surface.
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='jet',
                      linewidth=0, antialiased=False)
                  
# ------ Customize the z axis.
#ax.set_zlim(300, 400)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# ------- Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

# ------- Figure show-save
plt.show()
#plt.tight_layout()
#plt.savefig('/Users/ramytanios/Desktop/Computational Fluid Dynamics Codes/Fortran Codes/Diffusion Convection/results/addSmart.pdf')

