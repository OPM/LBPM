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

# The flag of writing data in binary format
# if True:  write data to binary format (besides the normal *.nc format)
# if False: not write data to binary format
#flag_write_bin = int(sys.argv[2])

# read Domain.in 
size_of_file=os.path.getsize('./00000') # the size of the semgented file in bytes
Size_of_Var = 8 #size_of(double) = 8 bytes
Num_of_Var = 4 # (Phase, Pressure, SignDist, BlobID)

#TODO: add input option so that user can choose what the output data they want, 
#      currently the output data is unnecessarily too big
#      perhaps just live phase, pressure and sign_dist
offset_data_header_1 = len("Var: domain-00000-Phase: 1, 3, "+str(lx*ly*lz)+", "+str(Size_of_Var*lx*ly*lz)+", double\n")
offset_data_header_2 = len("Var: domain-00000-Pressure: 1, 3, "+str(lx*ly*lz)+", "+str(Size_of_Var*lx*ly*lz)+", double\n")
offset_data_header_3 = len("Var: domain-00000-SignDist: 1, 3, "+str(lx*ly*lz)+", "+str(Size_of_Var*lx*ly*lz)+", double\n")
offset_data_header_4 = len("Var: domain-00000-BlobID: 1, 3, "+str(lx*ly*lz)+", "+str(Size_of_Var*lx*ly*lz)+", double\n")
offset_Var_with_endian = Num_of_Var*lx*ly*lz*Size_of_Var + Num_of_Var*1
offset_pre_header = size_of_file-(offset_Var_with_endian+offset_data_header_1+offset_data_header_2+offset_data_header_3+offset_data_header_4)

data_group = glob('00*') #NOTE Might need to change this if your original image is segmented into too many (e.g. hundreds of) parts
data_group.sort()
if not data_group:
    print 'Error: No data files: '+data_group
else:
    for ii in range(len(data_group)):
        print '**Info: Read data file: '+data_group[ii]
        phase=np.memmap(data_group[ii],dtype=np.float64,\
                        offset=offset_pre_header+offset_data_header_1,\
                        shape=lx*ly*lz)
        pressure=np.memmap(data_group[ii],dtype=np.float64,\
                           offset=offset_pre_header+offset_data_header_1+Size_of_Var*lx*ly*lz+len("\n")+offset_data_header_2,\
                           shape=lx*ly*lz)
        blob_id=np.memmap(data_group[ii],dtype=np.float64,\
                          offset=offset_pre_header+offset_data_header_1+offset_data_header_2+offset_data_header_3+offset_data_header_4+3*Size_of_Var*lx*ly*lz+3*len("\n"),\
                          shape=lx*ly*lz)
        
        print '**Info: Write data in normal binary format......'
        phase.tofile('Phase.'+data_group[ii])
        pressure.tofile('Pressure.'+data_group[ii])
        blob_id.tofile('BlobID.'+data_group[ii])
        print '**Info: Data extraction is completed.'
    #end for
#end if    















