#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from netCDF4 import Dataset

def write_NetCDF_file_py(file_name,data_name,data_IN,data_IN_dtype,lx,ly,lz):
	
	# Important !
	# Note the dataOUT should have the shape (lz,ly,lx)
	data_IN.shape=(lz,ly,lx) # make sure data_IN has the right shape
	# open a new netCDF file for writing.
	ncfile = Dataset(file_name,'w') 
	# create the output data.
	
	# create the x, y and z dimensions.
	ncfile.createDimension('x',lx)
	ncfile.createDimension('y',ly)
	ncfile.createDimension('z',lz)
	# create the variable (4 byte integer in this case)
	# first argument is name of variable, second is datatype, third is
	# a tuple with the names of dimensions.
	data = ncfile.createVariable(data_name,data_IN_dtype,('z','y','x'))
	data[:] = data_IN
        #data.voxel_unit = 'micrometer'
        #data.voxel_size=5.7
	ncfile.close()
        print '**Info: The *.nc file is written successfully !'
#end def

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

## -------- Setup output data information -------------##
# 1. For *.nc file
# 1.1 for ID.* files
id_data_output_file='ID_All.nc'
id_data_output_variable='segmented'
# 1.2 for SignDist.* files
dist_data_output_file='SignDist_All.nc'
dist_data_output_variable='signed_distance'
# 1.3 for Phase field files
phase_data_output_file='Phase_All.nc'
phase_data_output_variable='phase_field'
# 2. For *.bin files
id_data_output_file_bin='ID_All.bin'
dist_data_output_file_bin='SignDist_All.bin'
phase_data_output_file_bin='Phase_All.bin'

# The size of the subdomain
lx+=2 # consider the halo layer
ly+=2
lz+=2

id_group = glob('ID.*')
id_group.sort()
dist_group = glob('SignDist.*')
dist_group.sort()

# Reconstruct ID.xxxxx data
for z_axis in np.arange(nproc_z):
    for y_axis in np.arange(nproc_y):
        for x_axis in np.arange(0,nproc_x,2):
            if nproc_x%2==0: # nproc_x is an even number
                temp_x1 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                temp_x2 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.int8)
                temp_x1.shape = (lz,ly,lx)
                temp_x2.shape = (lz,ly,lx)
                temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
            else: # nproc_x is an odd number
                if nproc_x==1:
                    temp_x1 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                elif x_axis < nproc_x-2: # i.e. nproc_x = 3 or 5 or 7 ...
                    temp_x1 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                    temp_x2 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.int8)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x2.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                    temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                    temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                #end if
            #end if
            if x_axis == 0:
                ID_x_temp = temp_x1.copy()
            else:
                ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
            #end if
        #end for
        if nproc_x !=1 and nproc_x%2 != 0:
            temp_x1 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+nproc_x-1],dtype=np.int8)
            ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
        #end if

        if y_axis == 0:
            ID_y_temp = ID_x_temp.copy()
        else:
            ID_y_temp = np.concatenate((ID_y_temp,ID_x_temp),axis=1)
        #end if
    #end for

    if z_axis == 0:
        ID_z_temp = ID_y_temp.copy()
    else:
        ID_z_temp = np.concatenate((ID_z_temp,ID_y_temp),axis=0)
    #end if
#end for

## ---------------- ID decoding ------------------------- ##
# By convention, in LBPM-WIA, we have segmented image as:
# 0 -> solid phase
# 1 -> non-wetting phase
# 2 -> wetting phase
## ---------------- ID decoding ------------------------- ##

# By default, the inlet layer is non-wetting phase, and the outlet layer is wetting phase
# Also by default, the flow direction is along z-axis, and the inlet/outlet plane is in x-y plane
ID_z_temp[0,:] =1
ID_z_temp[-1,:]=2
ID_z_temp.tofile(id_data_output_file_bin)
write_NetCDF_file_py(id_data_output_file,id_data_output_variable,ID_z_temp,ID_z_temp.dtype.char,LX,LY,LZ)

# Create phase field based on the ID profiles
# By convention, in LBPM-WIA, we have segmented image as:
# 0 -> solid phase
# 1 -> non-wetting phase
# 2 -> wetting phase
phase = ID_z_temp.copy()
phase = phase.astype(np.float32)
#TODO: make this file read in Color.in file so phi_s can be assigned according to Color.in
phase[phase==0.0] = -0.55 # solid phase
phase[phase==1.0] = 1.0   # non-wetting phase
phase[phase==2.0] = -1.0  # wetting pahse
phase.tofile(phase_data_output_file_bin)
write_NetCDF_file_py(phase_data_output_file,phase_data_output_variable,phase,phase.dtype.char,LX,LY,LZ)

# Reconstruct SignDist.xxxxx data
for z_axis in np.arange(nproc_z):
    for y_axis in np.arange(nproc_y):
        for x_axis in np.arange(0,nproc_x,2):
            if nproc_x%2==0: # nproc_x is an even number
                temp_x1 = np.fromfile(dist_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                temp_x2 = np.fromfile(dist_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                temp_x1.shape = (lz,ly,lx)
                temp_x2.shape = (lz,ly,lx)
                temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
            else: # nproc_x is an odd number
                if nproc_x==1:
                    temp_x1 = np.fromfile(id_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                elif x_axis < nproc_x-2: # i.e. nproc_x = 3 or 5 or 7 ...
                    temp_x1 = np.fromfile(dist_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                    temp_x2 = np.fromfile(dist_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x2.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                    temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                    temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                #end if
            #end if
            if x_axis == 0:
                ID_x_temp = temp_x1.copy()
            else:
                ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
            #end if
        #end for
        if nproc_x !=1 and nproc_x%2 !=0:
            temp_x1 = np.fromfile(dist_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+nproc_x-1],dtype=np.float64)
            ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
        #end if

        if y_axis == 0:
            ID_y_temp = ID_x_temp.copy()
        else:
            ID_y_temp = np.concatenate((ID_y_temp,ID_x_temp),axis=1)
        #end if
    #end for

    if z_axis == 0:
        ID_z_temp = ID_y_temp.copy()
    else:
        ID_z_temp = np.concatenate((ID_z_temp,ID_y_temp),axis=0)
    #end if
#end for
ID_z_temp.tofile(dist_data_output_file_bin)
write_NetCDF_file_py(dist_data_output_file,dist_data_output_variable,ID_z_temp,ID_z_temp.dtype.char,LX,LY,LZ)


