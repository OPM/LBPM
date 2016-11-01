#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
#from glob import glob
import sys

# Check if there is a proper command line argument
#if len(sys.argv) !=2:
#    sys.stderr.write('Usage: ' + sys.argv[0] + ' Domain.in\n')
#    sys.exit()
## end if

# Read 'Domain.in' file
f = open('Domain.in','r') # read-only
lines = f.readlines()
nprocx, nprocy, nprocz = np.fromstring(lines[0].splitlines()[0],dtype=np.int32,sep=' ')
nx, ny, nz = np.fromstring(lines[1].splitlines()[0],dtype=np.int32,sep=' ')
Lx, Ly, Lz = np.fromstring(lines[3].splitlines()[0],dtype=np.float32,sep=' ')
f.close()


lx = 3
ly = 500
lz = 500

# Load the data
SignDist = np.fromfile('SignDist.00000',dtype = np.float64)
SignDist.shape = (lz+2,ly+2,lx+2)
ID = np.fromfile('ID.00000',dtype = np.uint8)
ID.shape = (lz+2,ly+2,lx+2)

# Plot
#plt.figure(1)
#plt.title('SignDist Map')
#plt.pcolormesh(SignDist[:,:,2])
#plt.grid(True)
#plt.axis('equal')
#plt.colorbar()
#
#
#plt.figure(2)
#plt.title('ID.xxxxx')
#plt.pcolormesh(ID[:,:,2],cmap='hot')
#plt.grid(True)
#plt.axis('equal')
#plt.show()





