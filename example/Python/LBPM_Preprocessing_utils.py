#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset


def write_NetCDF_file_py(file_name,data_name,data_IN,data_IN_dtype,lx,ly,lz):
	
    # NOTE: the dataIN should have the size of lz*ly*lx
	data_IN.shape=(lz,ly,lx) # make sure data_IN has the right size
	# open a new netCDF file for writing.
	ncfile = Dataset(file_name,'w') 
	# create the output data.
	
	# create the x, y and z dimensions.
	ncfile.createDimension('x',lx)
	ncfile.createDimension('y',ly)
	ncfile.createDimension('z',lz)

	# create the variable
	# first argument is name of variable, second is datatype, third is
	# a tuple with the names of dimensions.
	data = ncfile.createVariable(data_name,data_IN_dtype,('z','y','x'))
	data[:] = data_IN
#    data.voxel_unit = 'micrometer'
#    data.voxel_size=5.7
	ncfile.close()
        print '**Info: The *.nc file is written successfully !'
#end def

def read_NetCDF_file_py(file_name,data_name):

    ncfile = Dataset(file_name,mode='r')
    dataOUT = ncfile.variables[data_name][:]
    
    ncfile.close()
    return dataOUT
#end def



