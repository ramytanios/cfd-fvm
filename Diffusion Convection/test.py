from numpy import *
import matplotlib.pyplot as plt

it,er = genfromtxt(r'test.txt', unpack=True) 
plt.scatter(it, er)
plt.pause(0.05)

#plt.show()
