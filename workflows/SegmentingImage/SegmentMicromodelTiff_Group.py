#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from glob import glob
from PIL import Image # import Python Imaging Library (PIL)
import os
plt.rcParams["savefig.directory"] = os.getcwd() # change the default pylab 'savefig' path (i.e. ~/) to cwd
from SegmentMicromodelTiff_utils import convert_image_to_array,read_input_parameters

# Load a group of *.in files 
# e.g. 'parameters_1.in' etc.
file_input_group = glob('parameters_*.in')
#output_result = np.empty((0,),dtype=float)

if not file_input_group:
    print 'Error: Input files cannot be found ! '
else:
    for ii in range(len(file_input_group)):
        file_input_single = file_input_group[ii]
        Para = read_input_parameters(file_input_single)
        # set your output file names
        input_file_idx = file_input_single[file_input_single.find('parameters_')+len('parameters_'):file_input_single.find('.in')]
        output_file_name = 'Segmented_image_test_'+input_file_idx+'.dat'

        imin = Para['imin']
        imax = Para['imax']
        imin_b = Para['imin_b'] # greater than 'imin'
        imax_b = Para['imax_b'] # smaller than 'imax'
        top_layer = Para['top_layer']
        bot_layer = Para['bot_layer']
        DEPTH = Para['DEPTH']
        NY = imax - imin
        NZ = NY
        NX = DEPTH+top_layer+bot_layer

        print "Processing "+Para['image_file_name']+" now..."
        print "------ Micromodel dimensions (NX, NY, NZ) are (%i, %i, %i) ------" % (NX,NY,NZ)
        
        # parameters for segmentation
        threshold_nw = Para['threshold_nw'] # NW phase: RGB values > threshold_nw
        threshold_s  = Para['threshold_s']  # solid:    RGB values < threshold_s
        # W phase: threshold_s <= RGB values <= threshold_nw

        # Get an array from the input .tiff image
        im_array = convert_image_to_array(Para['image_file_name'],imin,imax)

        # Initialize the segmented image 'ID'
        # NOTE: 'ID' has the dimension (height, width, DEPTH)
        ID = np.zeros((NZ,NY,NX),dtype=np.uint8)
        # Map the picels to the 3D geometry
        # 1. Identify the non-wetting phase
        ID[im_array>threshold_nw,top_layer:top_layer+DEPTH] = 1
        # 2. Identify the wetting phase
        ID[np.logical_and(im_array>=threshold_s,im_array<=threshold_nw),\
           top_layer:top_layer+DEPTH] = 2
        # 3. Post boundary retreatment along y-axis. Note z-axis is always the flow direction
        ID[:,:imin_b-imin,top_layer:top_layer+DEPTH]=0
        ID[:,ID.shape[1]-(imax-imax_b):,top_layer:top_layer+DEPTH]=0

        # Calculate the porosity
        POROSITY = 1.0 - 1.0*ID[ID==0].size/ID.size
        print "Porosity of the micromodel is: "+str(POROSITY) 


        # Write out the segmented data (in binary format)
        print "The segmented image is written to "+output_file_name+" now..."
        ID.tofile(output_file_name)
        # NOTE when you want to load the binary file
        # Do as follows:
        # ID_reload = np.fromfile(file_name,dtype=np.uint8)

    #end for
#end if

    
# In case you want to test the segmentation
plt.figure(1)
#plt.subplot(1,2,1)
#plt.title('original RGB figure')    
#plt.pcolormesh(im_array);
#plt.axis('equal')
#plt.colorbar()
#plt.subplot(1,2,2)
plt.title('Segmented image')
# This will show the last-processed segmented image
cmap = plt.cm.get_cmap('hot',3) #Plot 3 discrete colors for NW, W and Solids
cax=plt.pcolormesh(ID[:,:,NX/2],cmap=cmap,vmin=-0.5,vmax=2.5);
cbar = plt.colorbar(cax,ticks=[0,1,2])
plt.axis('equal')
plt.show()




