import numpy as np
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
from glob import glob

def read_NetCDF_file_py(file_name,data_name):

    ncfile = Dataset(file_name,mode='r')
    dataOUT = ncfile.variables[data_name][:]
    
    ncfile.close()
    return dataOUT
#end def

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

#NOTE: To do a 3-phase segmentation, you need the original segmented image file, 
#      which contains 0->solid and 1->pore.
print "**Info: Need to read the original input segmetned image, "
print "        which contains 0 -> solid and 1 -> pore"
img_name = "domain1_64_reduced_lbpm.bin"
print "**Info: Read the image: "+img_name
img = np.fromfile(img_name,dtype = np.int8)

input_nc_file_group = glob('OutputData_vis*.nc')
input_nc_file_var_name          = 'Phase'
output_nc_file_var_name         = 'Phase_seg'


if not input_nc_file_group:
    print "**Error: No input data files: "+input_nc_file_group+" have been found !"
else:
    for ii in range(len(input_nc_file_group)):
        input_nc_file_name = input_nc_file_group[ii]
        output_nc_file_name = input_nc_file_name[:-3]+'_PhaseSeg.nc'
        print "**Info: Read data file: "+input_nc_file_name
        phase=read_NetCDF_file_py(input_nc_file_name,input_nc_file_var_name)
        (lz,ly,lx) = phase.shape
        img.shape = phase.shape

        # Double check the phase field has the same size as the segmented image
        if phase.size != img.size:
            print "**Error: The size of the loaded segmented image does not match the size of the phase field data !"
        else:
            phase[phase>0.0]=1
            phase[phase<0.0]=2
            phase[img==0]=0
            phase = phase.astype(np.int32)
            write_NetCDF_file_py(output_nc_file_name,output_nc_file_var_name,phase,phase.dtype.char,lx,ly,lz)
            print "**Info: The segmented phase field from the floating phase field has been written."
        #end if
    #end for
#end if    








