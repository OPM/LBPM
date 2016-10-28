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

#def get_segmented_array(image_array,parameters,NX,NY,NZ):
#    # 'image_array' is a numpy array of RGB values, with the same size as micromodel 
#    # 'parameters' is a dictionary of input parameters
#    
#    # Initialize the segmented image 'ID'
#    # NOTE: 'ID' has the dimension (DEPTH, height, width)
#    ID = np.zeros((NX,NZ,NY),dtype=np.uint8)
#    # Map the picels to the 3D geometry
#    # 1. Identify the non-wetting phase
#    ID[top_layer:top_layer+DEPTH,im_array>threshold_nw] = 1
#    # 2. Identify the wetting phase
#    ID[top_layer:top_layer+DEPTH,\
#    np.logical_and(im_array>=threshold_s,im_array<=threshold_nw)] = 2
#    # 3. Post boundary retreatment along y-axis
#    ID[top_layer:top_layer+DEPTH,:,:imin_b-imin]=0
#    ID[top_layer:top_layer+DEPTH,:,ID.shape[2]-(imax-imax_b):]=0
#    return ID
##end def
