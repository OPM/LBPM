#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from glob import glob
from PIL import Image # import Python Imaging Library (PIL)
import sys
from SegmentMicromodelTiff_utils import *

# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('Usage: ' + sys.argv[0] + ' <Input parameter file>\n')
    print 'Input parameter file: '
    print 'line 0: the base name of the experimental images\n',\
          'line 1: imin imax\n',\
          'line 2: imin_b imax_b\n',\
          'line 3: top_layer bottom_layer DEPTH\n',\
          'line 4: threshold_nw threshold_s'
    sys.exit()
# end if

input_parameter_file = sys.argv[1] # Give the name of the input parameter file 
BaseName_output_init_config = '_InitConfig_segmented.raw'# The string will be appended at the end of the 
                                                         # name of the input image as an output file name
BaseName_output_domain = '_segmented.raw'# The string will be appended at the end of the 
                                         # name of the input image as an output file name
image_format = '.tiff' # format of experimental images

# Read other input parameters to get NX, NY and NZ
Para = read_input_parameters(input_parameter_file)
NY = Para['imax'] - Para['imin']
NZ = NY
NX = Para['DEPTH']+Para['top_layer']+Para['bot_layer']

# Load a group of image files with 'base_name' in the names of files 
# e.g. 'Micromodel_1.tiff' etc.
file_input_group = glob('*'+Para['base_name']+'*'+image_format) # need to match the image format

# Process all imported experimental images
if not file_input_group:
    print 'Error: Input files cannot be found ! '
else:
    for ii in range(len(file_input_group)):
        file_input_single = file_input_group[ii]
        print "Processing image "+file_input_single+" now..."
        print "------ Micromodel dimensions (NX, NY, NZ) are (%i, %i, %i) ------" % (NX,NY,NZ)
        
        # Get an array from the input .tiff image
        im_array = convert_image_to_array(file_input_single,Para['imin'],Para['imax'])
        # Get initial fluid configuration shown in the *.Tiff image
        ID = get_initial_config(im_array,Para,NX,NY,NZ)
        # Get only the domain (0: solid, 1: fluid) 
        PorousMedium = get_porous_medium(im_array,Para,NX,NY,NZ)
        # Calculate the porosity
        POROSITY = 1.0 - 1.0*ID[ID==0].size/ID.size
        print "Porosity of the micromodel is: "+str(POROSITY) 
        
        # Write the segmeted data to the output
        # set the output file names (need to chop off say the '.tiff' part)
        output_file_name_domain      = file_input_single[:file_input_single.find('.')]+BaseName_output_domain
        output_file_name_init_config = file_input_single[:file_input_single.find('.')]+BaseName_output_init_config
        # Write out the segmented data (in binary format)
        print "The segmented fluid configuration is written to "+output_file_name_init_config+" now..."
        print "The segmented porous domain is written to "+output_file_name_domain+" now..."
        ID.tofile(output_file_name_init_config)
        PorousMedium.tofile(output_file_name_domain)
        # NOTE when you want to load the binary file
        # Do as follows:
        # ID_reload = np.fromfile(file_name,dtype=np.uint8)

    #end for
#end if

    
### In case you want to test the segmentation
#plt.figure(1)
#plt.subplot(1,2,1)
#plt.title('original RGB figure')    
#plt.pcolormesh(im_array);
#plt.axis('equal')
#plt.colorbar()
#plt.subplot(1,2,2)
#plt.title('Segmented image')
### This will show the last-processed segmented image
#cmap = plt.cm.get_cmap('hot',3) #Plot 3 discrete colors for NW, W and Solids
#cax=plt.pcolormesh(ID[:,:,NX/2],cmap=cmap,vmin=-0.5,vmax=2.5);
#cbar = plt.colorbar(cax,ticks=[0,1,2])
#plt.axis('equal')
#plt.show()




