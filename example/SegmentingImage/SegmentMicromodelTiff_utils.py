#!/usr/bin/env python

import numpy as np
from PIL import Image # import Python Imaging Library (PIL)

def convert_image_to_array(file_name,imin,imax):
    # file_name: the name of the input image (a string)
    # imin: (left) boundary that defines the output array
    # imax: (right) boundary that defines the output array
    # Note: assume a grayscale image (whose R, G, and B values are the same)
    # Return an numpy array with RGB values that has the same size as the wanted micromodel
    im = Image.open(file_name).convert('RGB')
    im_array = np.array(im) # 'im_array' has dimension (height, width, 3)
    # Grayscale input image is assumed, thus its R, G, and B values are the same
    im_array = im_array[:,:,0] # 2D array with dimension (height, width)
    im_array = im_array[:,imin:imax] # Chop the original image to a desirable size
    #im_array = np.flipud(im_array[:,imin:imax]) # You might want to flip the Z-axis
    return im_array
#end def

def read_input_parameters(input_file_name):
    # input_file_name is a string
    # The *.in file has the following lines
    # line 0: the base name of the experimental images
    # line 1: imin imax
    # line 2: imin_b imax_b
    # line 3: top_layer bottom_layer DEPTH
    # line 4: threshold_nw threshold_s
    # Note:  'threshold_nw' means: NW phase: RGB values > threshold_nw
    #        'threshold_s'  means: solid:    RGB values < threshold_s
    f = open(input_file_name,'r') # read-only
    lines = f.readlines()
    output_file = {'base_name':lines[0].splitlines()[0]}
    line1_array = np.fromstring(lines[1].splitlines()[0],dtype=np.int32,sep=' ')
    line2_array = np.fromstring(lines[2].splitlines()[0],dtype=np.int32,sep=' ')
    line3_array = np.fromstring(lines[3].splitlines()[0],dtype=np.int32,sep=' ')
    line4_array = np.fromstring(lines[4].splitlines()[0],dtype=np.int32,sep=' ')
    output_file['imin']=line1_array[0]
    output_file['imax']=line1_array[1]
    output_file['imin_b']=line2_array[0]
    output_file['imax_b']=line2_array[1]
    output_file['top_layer']=line3_array[0]
    output_file['bot_layer']=line3_array[1]
    output_file['DEPTH']=line3_array[2]
    output_file['threshold_nw']=line4_array[0]
    output_file['threshold_s'] =line4_array[1]
    f.close()
    return output_file
#end def

def get_initial_config(image_array,Para,NX,NY,NZ):
    # This function will give the same fluid configuration shown in the *.Tiff image
    # Function input:
    # 'image_array' is a numpy array of RGB values, with the same size as micromodel 
    # 'Para'        is a dictionary of input parameters

    # Initialize the segmented image 'ID'
    # NOTE: 'ID' has the dimension (height, width, DEPTH) or (NZ,NY,NX)
    ID = np.zeros((NZ,NY,NX),dtype=np.uint8)
    # Map the pixels to the 3D geometry
    # Rules: NW phase: RGB values > threshold_nw
    #        Solid:    RGB values < threshold_s
    # 1. Identify the non-wetting phase
    ID[image_array>Para['threshold_nw'],Para['top_layer']:Para['top_layer']+Para['DEPTH']] = 1
    # 2. Identify the wetting phase
    ID[np.logical_and(image_array>=Para['threshold_s'],image_array<=Para['threshold_nw']),\
       Para['top_layer']:Para['top_layer']+Para['DEPTH']] = 2
    # 3. Post boundary retreatment along y-axis. Note z-axis is always the flow direction
    ID[:,:Para['imin_b']-Para['imin'],Para['top_layer']:Para['top_layer']+Para['DEPTH']]=0
    ID[:,ID.shape[1]-(Para['imax']-Para['imax_b']):,Para['top_layer']:Para['top_layer']+Para['DEPTH']]=0
    return ID
#end def    

def get_porous_medium(image_array,Para,NX,NY,NZ):
    # This function only gives the domain (0->solid,1->fluid) of the *.Tiff image
    # Function input:
    # 'image_array' is a numpy array of RGB values, with the same size as micromodel 
    # 'Para'        is a dictionary of input parameters

    domain = np.zeros((NZ,NY,NX),dtype=np.uint8)
    # Map the pixels to the 3D geometry
    # Rules: NW phase: RGB values > threshold_nw
    #        Solid:    RGB values < threshold_s
    # 1. Identify fluid phase
    domain[image_array>=Para['threshold_s'],Para['top_layer']:Para['top_layer']+Para['DEPTH']]=1
    # 2. Post boundary retreatment along y-axis. Note z-axis is always the flow direction
    domain[:,:Para['imin_b']-Para['imin'],Para['top_layer']:Para['top_layer']+Para['DEPTH']]=0
    domain[:,domain.shape[1]-(Para['imax']-Para['imax_b']):,Para['top_layer']:Para['top_layer']+Para['DEPTH']]=0
    return domain
#end def

