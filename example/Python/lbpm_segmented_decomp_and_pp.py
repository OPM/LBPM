#!/usr/bin/env python
import sys
import numpy as np
import skfmm
#NOTE: need python package 'array_split'
#      pip install array-split
#      Ref: https://pypi.python.org/pypi/array-split/0.1.3
#      For more information on using array_split function, please see
#      http://array-split.readthedocs.io/en/latest/reference/generated/array_split.array_split.html#array_split.array_split
from array_split import array_split 

# Check if there is a proper command line argument
if len(sys.argv) !=3:
    sys.stderr.write('**Error: Usage: ' + sys.argv[0] + '<Segmented.in> <Domain.in>\n')
    sys.exit()
# end if

## ***************** Read 'Segmented.in' ******************** ##
print "**Info: Reading Segmented.in file."
f = open(sys.argv[1],'r')
lines = f.readlines()
# Read the file name of the segmented image
seg_data_name=lines[0]
seg_data_name=seg_data_name[:-1]
# Read the (original) size of the segmeted image
LX,LY,LZ = np.fromstring(lines[1].splitlines()[0],dtype = np.int32,sep=' ') 
f.close()

## ******************************* Read 'Domain.in' ************************************* ##
print "**Info: Reading Domain.in file."
f = open(sys.argv[2],'r')
lines = f.readlines()
# Read processors configuration
nproc_x,nproc_y,nproc_z=np.fromstring(lines[0].splitlines()[0],dtype = np.int32,sep=' ')
# Read subdomain dimensions
lx,ly,lz = np.fromstring(lines[1].splitlines()[0],dtype = np.int32,sep=' ')
# Read domain dimensions
LX_temp,LY_temp,LZ_temp = np.fromstring(lines[3].splitlines()[0],dtype = np.float32,sep=' ') 
LX_temp=int(LX_temp)
LY_temp=int(LY_temp)
LZ_temp=int(LZ_temp)
f.close()
# Double check if the input of 'Domain.in' is consistent with that in the 'Segmented.in'    
if LX != LX_temp or LY != LY_temp or LZ != LZ_temp:
    print "**Error: Inconsistent image size input between Domain.in and Segmented.in !"
    print "         The following domain decomposition may fail due to this inconsistency."
#end if

print "**Info: Start preparing the domain decomposition and ID field labelling..."
print "**Info: The node type convention for LBPM-WIA is as follows:"
print "*********************************"
print "**    0 -> solid phase         **"
print "**    1 -> non-wetting phase   **"
print "**    2 -> wetting phase       **"
print "*********************************"

# Load the segmented image
# NOTE: Assume the input segmented image has ID field as follows:
# 0 -> solid nodes
# 1 or 2 -> fluid nodes
seg_data = np.fromfile(seg_data_name,dtype=np.int8)
#TODO: Make a function for calculating the signed distance, so to save the memory
seg_data_domain = seg_data.copy() # seg_data_domain is used for calculating the signed distance
seg_data_domain[seg_data_domain == 0] = -1
seg_data_domain[seg_data_domain  > 0] =  1
seg_data.shape = (LZ,LY,LX)
seg_data_dist = skfmm.distance(seg_data_domain.reshape((LZ,LY,LX)))

# Decomposition
seg_data_decomp_list      = array_split(seg_data,axis=[nproc_z,nproc_y,nproc_x])
seg_data_dist_decomp_list = array_split(seg_data_dist,axis=[nproc_z,nproc_y,nproc_x])
#TODO: check the size of the decomposed data, comparing to the loaded subdomain size (lz,ly,lx)

# Write the decomposed data
ID_prefix="ID."
SignDist_prefix="SignDist."
halo_layer = 1 #The length of the halo layer
for ii in np.arange(nproc_z*nproc_y*nproc_x):
    ID_name=ID_prefix+"%05i" % ii
    SignDist_name=SignDist_prefix+"%05i" % ii
    temp = seg_data_decomp_list[ii]
    temp = np.lib.pad(temp,((halo_layer,halo_layer),(halo_layer,halo_layer),(halo_layer,halo_layer)),'edge')
    print "**Info: Write to file: "+ID_name
    temp.tofile(ID_name)
    temp = seg_data_dist_decomp_list[ii]
    temp = np.lib.pad(temp,((halo_layer,halo_layer),(halo_layer,halo_layer),(halo_layer,halo_layer)),'edge')
    temp.tofile(SignDist_name)
    print "**Info: Write to file: "+SignDist_name
#end for







