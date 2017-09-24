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
# -------------------------------------------------------------------------------------------

phase=Reconstruct_single_data('Phase.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))
pressure=Reconstruct_single_data('Pressure.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))
blob_id=Reconstruct_single_data('BlobID.0*',(nproc_x,nproc_y,nproc_z),(lx,ly,lz))
output_data=np.hstack((phase,pressure,blob_id))

#output_data_name_bin='Output_data.bin'
folder_path=os.getcwd()
output_data_name_nc ='OutputData_'+os.path.basename(folder_path)+'.nc'
# Write to binary file
#output_data.tofile(output_data_name_bin)
#print '**Info: Output data \''+output_data_name_bin+'\' is written.'
# Write to *nc file
Writer_NetCDF_multiVar_LBPMWIA(output_data_name_nc,output_data.astype(np.float32),LX,LY,LZ)

# Clean up temporary files
files=glob('Phase.0*')
for file_tmp in files:
    print "**Warning: Remove file: "+file_tmp
    os.remove(file_tmp)
#end for
files=glob('Pressure.0*')
for file_tmp in files:
    print "**Warning: Remove file: "+file_tmp
    os.remove(file_tmp)
#end for
files=glob('BlobID.0*')
for file_tmp in files:
    print "**Warning: Remove file: "+file_tmp
    os.remove(file_tmp)
#end for



