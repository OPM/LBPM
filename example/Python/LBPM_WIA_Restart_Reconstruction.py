#!/usr/bin/env python

import sys
import os
import numpy as np
#import matplotlib.pyplot as plt
from glob import glob
#from netCDF4 import Dataset
from LBPM_WIA_OutputData_Postprocessing_utils import *


# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('**Error: Usage: ' + sys.argv[0] + ' <Domain.in>\n')
    sys.exit()
# end if

# Read 'Domain.in' --------------------------------------------------------------------------
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

# NOTE: Do not forget the halo layer attached on each subdomain !
lx = lx + 2
ly = ly + 2
lz = lz + 2
# -------------------------------------------------------------------------------------------

print "**Info: Reconstructing individual data sets according to the processor grid..."
dataOUT = Reconstruct_SignDist_data('DenNW.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('DenW.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f0.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f1.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f2.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f3.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f4.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f5.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f6.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f7.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f8.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f9.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f10.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f11.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f12.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f13.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f14.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f15.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f16.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f17.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
dataOUT = np.hstack((dataOUT,Reconstruct_SignDist_data('f18.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))))
#output_data_name_bin='Output_data.bin'
output_data_name_nc ='Restart_dbg_TimeStep_.nc'
# Write to binary file
#output_data.tofile(output_data_name_bin)
#print '**Info: Output data \''+output_data_name_bin+'\' is written.'
# Write to *nc file
Writer_NetCDF_multiVar_Restart_LBPMWIA(output_data_name_nc,dataOUT.astype(np.float32),LX,LY,LZ)

# Clean up temporary files
files=glob('Den*')
for file_tmp in files:
    print "**Warning: Remove file: "+file_tmp
    os.remove(file_tmp)
#end for
files=glob('f*.0*')
for file_tmp in files:
    print "**Warning: Remove file: "+file_tmp
    os.remove(file_tmp)
#end for



