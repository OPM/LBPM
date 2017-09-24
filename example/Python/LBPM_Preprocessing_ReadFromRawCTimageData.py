#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from LBPM_Preprocessing_utils import *
import matplotlib.pyplot as plt

# Load the *nc image
# NOTE: for the data name, you might want to check the *.nc file with ncdump -h
image = read_NetCDF_file_py('subdomain_1_from_Anna.nc','segmented')
# The loaded data from the original CT image is a masked array
# with the pore space being masked.
image = np.ma.filled(image,fill_value=0)
(lz,ly,lx) = image.shape
# NOTE: the CT image has the convention as:
#       0 -> pore space
#       1 -> solid phase
# Whereas LBPM-WIA has the convention as:
# 0 -> solid phase
# 1 -> non-wetting phase
# 2 -> wetting phase
image = np.logical_not(image)
image = image.astype(np.int8)# Now: 0 -> solid phase
                             #      1 -> fluid phase
image_red = image[0:256,0:256,0:256]
(lz_red,ly_red,lx_red) = image_red.shape

## 0.0. Save the original chopped image, in simple binary format.
#image_red.tofile('subdomain_1.bin')
## 0.1. Save the original chopped image, in simple *nc format
write_NetCDF_file_py('subdomain_1.nc','segmented',\
                     image_red,image_red.dtype.char,lx_red,ly_red,lz_red)
image_red.tofile('subdomain_1.bin')
# Preprocess the image
# ---------------- By default: ---------------------------------------------
# 1. The flow direction is along positive z-axis
#    z_min: inlet;  z_max: outlet
# 2. The boundary condition, e.g. pressure boundary condition is applied
#    on x-y planes
# --------------------------------------------------------------------------

# Steps of preprocessing:
# 1. wrap boundaries along z-axis with solid nodes
# 1.1 (optional) save the raw image for morphological drainge & opeing
# 2. change the original z_min and z_max layers to reservoir layers
# 3. set z_min reservoir layer to non-wetting phase reservoir
# 4. set z_max reservoir layer to wetting phase reservoir
# 5. add porous plate layer before the wetting phase reservoir
# 6. saturate the pore space with wetting phase

# 1. wrap the boundaries
# TODO: make a variable for wall_width
image_red[:,0,:]        = 0
image_red[:,ly_red-1,:] = 0
image_red[:,:,0]        = 0
image_red[:,:,lx_red-1] = 0

# 1.1 (optional) Save the image for morphological analysis
#print "**Info: Save data for morphological drainage and imbibition."
#image_red.tofile('subdomain_1_reduced_morph.bin') # dont use this

# 2. Saturate the medium with wetting phase
image_red[image_red==1] = 2
# 3. Set up reservoir layers
# Input parameters for the reservoir layers
NWR_len = 1 # The thickness of the non-wetting phaer reservoir
WR_len  = 1 # The thickness of the wetting phase reservoir
image_red[0:NWR_len,:]  = 1
image_red[-WR_len:,:]   = 2
# temp test
#plt.pcolormesh(image_red[:,1,:],cmap='hot')
#plt.axis('equal')
#plt.colorbar();
#plt.show()
# 4. setup up a porous plate saturated with **wetting phase**
# Input parameters for the porous plate
pp_size = 3 #NOTE: Although it is made as an variable, pp_size of 3 is recommended
pp_len  = 3 # the thickness of the porous plate
pp = 2*np.ones((pp_size,pp_size),dtype=np.int8)
#pp = np.lib.pad(pp,((1,1),(1,1)),'constant',constant_values=0)
pp = np.lib.pad(pp,((1,1),(1,1)),'constant',constant_values=0)
(ly_pp,lx_pp) = pp.shape
pp_pad = np.tile(pp,(ly_red/ly_pp,lx_red/lx_pp))
l_residual=ly_red%ly_pp
pp_pad = np.lib.pad(pp_pad,((l_residual%2,0),(l_residual%2,0)),'constant',constant_values=0)
pp_pad = np.lib.pad(pp_pad,((l_residual/2,l_residual/2),(l_residual/2,l_residual/2)),'constant',constant_values=0)
# 5. Create the simulation domain
image_red[-(WR_len+pp_len):-WR_len,:,:]=pp_pad
# temp test
#plt.pcolormesh(image_red[-5,:],cmap='hot')
#plt.axis('equal')
#plt.grid(True)
#plt.colorbar()
#plt.show()
# write to *.nc file for futher check
#print "**Info: Save data for LBPM simulations."
#write_NetCDF_file_py('subdomain_1_reduced_lbpm.nc','segmented',\
#                     image_red,image_red.dtype.char,lx_red,ly_red,lz_red)
#image_red.tofile('subdomain_1_reduced_lbpm.bin')




