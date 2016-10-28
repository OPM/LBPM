#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt

lx = 3
ly = 500
lz = 500

# Load the data
SignDist = np.fromfile('SignDist.00000',dtype = np.float64)
SignDist.shape = (lz+2,ly+2,lx+2)
ID = np.fromfile('ID.00000',dtype = np.uint8)
ID.shape = (lz+2,ly+2,lx+2)

# Plot
plt.figure(1)
plt.title('SignDist Map')
plt.pcolormesh(SignDist[:,:,2])
plt.grid(True)
plt.axis('equal')
plt.colorbar()


plt.figure(2)
plt.title('ID.xxxxx')
plt.pcolormesh(ID[:,:,2],cmap='hot')
plt.grid(True)
plt.axis('equal')
plt.show()





