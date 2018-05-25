#!/usr/bin/env python

import sys
import os
import numpy as np
from glob import glob
#from netCDF4 import Dataset
#import matplotlib.pyplot as plt

# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('**Error: Usage: ' + sys.argv[0] + ' <Domain.in>\n')
    sys.exit()
# end if

# Read 'Domain.in' to obtain the size of 'ID.xxxxxx' and 'SignDist.xxxxx'
f = open(sys.argv[1],'r')
lines = f.readlines()
# Read processors configuration
nproc_x,nproc_y,nproc_z=np.fromstring(lines[0].splitlines()[0],dtype = np.int32,sep=' ')
# Read subdomain dimensions
lx,ly,lz = np.fromstring(lines[1].splitlines()[0],dtype = np.int32,sep=' ')
# Read domain dimensions
LX,LY,LZ = np.fromstring(lines[3].splitlines()[0],dtype = np.float32,sep=' ') 
LX=int(LX)
LY=int(LY)
LZ=int(LZ)
f.close()

# The total size of the subdomain including the halo layer
N=(lx+2)*(ly+2)*(lz+2)

# Read 'Restart.00*' binary files
# NOTE: Currently there are in total 21 sets of data in the Restart files
# They are: { rho_nw, rho_w, f_Even, f_odd }
# f_even -> { f0, f2, f4, f6, f8, f10, f12, f14, f16, f18 }
# f_odd  -> { f1, f3, f5, f7, f9, f11, f13, f15, f17 }
data_group = glob('Restart.00*')
data_group.sort()
if not data_group:
    print 'Error: No data files: '+data_group
else:
    for ii in range(len(data_group)):
        print '**Info: Read data file: '+data_group[ii]
        Restart=np.fromfile(data_group[ii],dtype=np.float64)
        Restart.shape = (N,21)
        print '**Info: Write data in normal binary format......'
        Restart[:,0].tofile('DenNW.'+data_group[ii][len('Restart.'):])
        Restart[:,1].tofile('DenW.'+data_group[ii][len('Restart.'):])
        Restart[:,2].tofile('f0.'+data_group[ii][len('Restart.'):])
        Restart[:,3].tofile('f2.'+data_group[ii][len('Restart.'):])
        Restart[:,4].tofile('f4.'+data_group[ii][len('Restart.'):])
        Restart[:,5].tofile('f6.'+data_group[ii][len('Restart.'):])
        Restart[:,6].tofile('f8.'+data_group[ii][len('Restart.'):])
        Restart[:,7].tofile('f10.'+data_group[ii][len('Restart.'):])
        Restart[:,8].tofile('f12.'+data_group[ii][len('Restart.'):])
        Restart[:,9].tofile('f14.'+data_group[ii][len('Restart.'):])
        Restart[:,10].tofile('f16.'+data_group[ii][len('Restart.'):])
        Restart[:,11].tofile('f18.'+data_group[ii][len('Restart.'):])
        Restart[:,12].tofile('f1.'+data_group[ii][len('Restart.'):])
        Restart[:,13].tofile('f3.'+data_group[ii][len('Restart.'):])
        Restart[:,14].tofile('f5.'+data_group[ii][len('Restart.'):])
        Restart[:,15].tofile('f7.'+data_group[ii][len('Restart.'):])
        Restart[:,16].tofile('f9.'+data_group[ii][len('Restart.'):])
        Restart[:,17].tofile('f11.'+data_group[ii][len('Restart.'):])
        Restart[:,18].tofile('f13.'+data_group[ii][len('Restart.'):])
        Restart[:,19].tofile('f15.'+data_group[ii][len('Restart.'):])
        Restart[:,20].tofile('f17.'+data_group[ii][len('Restart.'):])
        print '**Info: Data extraction is completed.'
    #end for
#end if    


