#!/usr/bin/env python
# TODO: double check the issue of the "view of the array" for all the following functions
#       make sure it is the copy of the array, not the view of the array is used for any
#       subprocesses.
import numpy as np
import scipy.stats as stats
import scipy.ndimage.morphology as morphology
import scipy.ndimage.measurements as measurements
#from dist_func_utils import *

def generate_morph_imbib_config(seg_image_input,R_critical):
    # This funciton works well at high-Pc-low-Sw 
    # There is no retarded Sw when the curvature is high
    # The method of this functions:
    # 1. Perform erosion with the radius of R_critical
    # 2. Perform dilation with the radius of R_critical
    if seg_image_input.ndim ==2 :
        seg_image = seg_image_input>0.0
        pore_vol = 1.0*seg_image.sum()
        radius = R_critical

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)
        
        # Step 1.2: Perform erosion with radius of R_critical
        seg_image_ero = morphology.binary_erosion(seg_image,structure=circle,border_value=1) 
        # NOTE: Be careful with the 'border_value' of the erosion - should be 'border_value=1'

        # Step 2: Perform dilation with radius of R_critical
        seg_image_dil = morphology.binary_dilation(seg_image_ero,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of the array after dilation is 'bool'
        seg_image_dil = seg_image_dil.astype(np.uint8)
        seg_image_dil[np.logical_not(seg_image_dil)]=2
        seg_image_dil[np.logical_not(seg_image)] = 0

        Sw = (seg_image_dil==2).sum()/pore_vol
    else: # 3D
        seg_image = seg_image_input>0.0
        pore_vol = 1.0*seg_image.sum()
        radius = R_critical

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 + (grid[2]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)
        
        # Step 1.2: Perform erosion with radius of R_critical
        seg_image_ero = morphology.binary_erosion(seg_image,structure=circle,border_value=1) 
        # NOTE: Be careful with the 'border_value' of the erosion - should be 'border_value=1'

        # Step 2: Perform dilation with radius of R_critical
        seg_image_dil = morphology.binary_dilation(seg_image_ero,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of the array after dilation is 'bool'
        seg_image_dil = seg_image_dil.astype(np.uint8)
        seg_image_dil[np.logical_not(seg_image_dil)]=2
        seg_image_dil[np.logical_not(seg_image)] = 0

        Sw = (seg_image_dil==2).sum()/pore_vol
    #end if
    return (seg_image_dil,Sw)
#end def    

def generate_morph_imbib_curv(seg_image_input,R_critical):
    # This funciton works well at high-Pc-low-Sw 
    # There is no retarded Sw when the curvature is high
    # The method of this functions:
    # 1. Perform erosion with the radius of R_critical
    # 2. Perform dilation with the radius of R_critical
    # *********************************************************************
    # Input: seg_image: a well shaped segmented image with size (lz,ly,lx)
    #        seg_image has values as : NW phase -> 1
    #                                   W phase -> 2
    #                               solid phase -> 0
    # *********************************************************************
    if seg_image_input.ndim == 2:
        pore_vol = 1.0*(seg_image_input>0.0).sum()
        radius = R_critical
        print 'Morphological Imbibition: processing critical radius: '+str(radius)+' now......'

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)
        
        # Step 1.2: Perform erosion with radius of R_critical
        seg_image_ero = morphology.binary_erosion(seg_image_input>0.0,structure=circle,border_value=1) 
        # NOTE: Be careful with the 'border_value' of the erosion - should be 'border_value=1'

        # Step 2: Perform dilation with radius of R_critical
        seg_image_ero_dil = morphology.binary_dilation(seg_image_ero,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of the array after dilation is 'bool'
        seg_image_ero_dil[seg_image_input<=0.0]=False

        Sw = 1.0 - seg_image_ero_dil.sum()/pore_vol
    else:
        pore_vol = 1.0*(seg_image_input>0.0).sum()
        radius = R_critical
        print 'Morphological Imbibition: processing critical radius: '+str(radius)+' now......'

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 + (grid[2]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)
        
        # Step 1.2: Perform erosion with radius of R_critical
        seg_image_ero = morphology.binary_erosion(seg_image_input>0.0,structure=circle,border_value=1) 
        # NOTE: Be careful with the 'border_value' of the erosion - should be 'border_value=1'

        # Step 2: Perform dilation with radius of R_critical
        seg_image_ero_dil = morphology.binary_dilation(seg_image_ero,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of the array after dilation is 'bool'
        seg_image_ero_dil[seg_image_input<=0.0]=False

        Sw = 1.0 - seg_image_ero_dil.sum()/pore_vol
    #end if
    return Sw
#end def    

def generate_morph_drain_config(seg_image_input,R_critical):
    # The method for this function follows Hilper & Miller AWR(2001)
    # 1. Perform erosion for the pore space with radius of R_critical
    # 2. Label the eroded pore space, and leave only the pore space that is still 
    #    connected with the non-wetting phase reservoir
    # 3. Perform the dilation for the labelled pore space with radius of R_critical
    # **************************************************************************
    # Input: seg_image: a well shaped segmented image with size (lz,ly,lx)
    #        seg_image has values as : NW phase -> 1
    #                                   W phase -> 2
    #                               solid phase -> 0
    # **************************************************************************
    if seg_image_input.ndim ==2:

        seg_image = seg_image_input>0.0
        pore_vol = 1.0*seg_image.sum()
        radius = R_critical

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)

        # Step 1.2: Perform erosion on the pore space
        # NOTE: the dtype of 'seg_im_ero' is 'bool'
        seg_im_ero = morphology.binary_erosion(seg_image,structure=circle,border_value=1)
        # NOTE: 'border_value' for erosion should be 'True'

        # Step 2: Label the eroded pore space
        # NOTE: Assume the NW phase reservoir is at the first layer of the domain
        #       i.e. at seg_image[0,:] - adjust it if this does not suit your need
        # For erosion, assume that diagonals are not considered
        # For erosion, assume that diagonals are not considered
        seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(2,1))
        #seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(2,2))
        # NOTE: Here I assume the inlet is at the first layer of the array's axis=2 (i.e. domain[0,:,:])\
        #       You can always change to any other layers as the inlet for this drainage.
        label_check = seg_im_ero_label_temp[0,seg_im_ero_label_temp[0,:]!=0]
        label_check = np.unique(label_check)

        # NOTE the following lines are only for your to check things
        # ******************** For check *******************************#
        # It assign the labelled array as: NW -> 1, W -> 2, Solid -> 0
        #seg_im_ero_label_show = seg_im_ero_label.copy()
        #seg_im_ero_label_show[seg_im_ero_label_show !=1] = 2
        #seg_im_ero_label_show[np.logical_not(seg_image_2d)]=0
        # ******************** End: for check **************************#
        
        seg_im_ero_label = np.zeros_like(seg_im_ero_label_temp,dtype=bool)
        for labels in label_check:
            seg_im_ero_label = np.logical_or(seg_im_ero_label,seg_im_ero_label_temp==labels)
        seg_im_ero_label = seg_im_ero_label.astype(np.uint8)
        
        # Step 3: perform dilation on the labelled pore space 
        seg_im_ero_label_dil = morphology.binary_dilation(seg_im_ero_label,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of 'seg_im_ero_label_dil' is 'bool'
        seg_im_ero_label_dil = seg_im_ero_label_dil.astype(np.uint8)
        seg_im_ero_label_dil[np.logical_not(seg_im_ero_label_dil)]=2
        seg_im_ero_label_dil[np.logical_not(seg_image)]=0

        Sw = (seg_im_ero_label_dil==2).sum()/pore_vol
    else: # 3D porous medium
        
        seg_image = seg_image_input>0.0
        pore_vol = 1.0*seg_image.sum()
        radius = R_critical

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 + (grid[2]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)

        # Step 1.2: Perform erosion on the pore space
        # NOTE: the dtype of 'seg_im_ero' is 'bool'
        seg_im_ero = morphology.binary_erosion(seg_image,structure=circle,border_value=1)
        # NOTE: 'border_value' for erosion should be 'True'

        # Step 2: Label the eroded pore space
        # NOTE: Assume the NW phase reservoir is at the first layer of the domain
        #       i.e. at seg_image[0,:] - adjust it if this does not suit your need
        # For erosion, assume that diagonals are not considered
        seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(3,1))
        #seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(3,3))
        # NOTE: Here I assume the inlet is at the first layer of the array's axis=2 (i.e. domain[0,:,:])\
        #       You can always change to any other layers as the inlet for this drainage.
        label_check = seg_im_ero_label_temp[0,seg_im_ero_label_temp[0,:]!=0]
        label_check = np.unique(label_check)

        # NOTE the following lines are only for your to check things
        # ******************** For check *******************************#
        # It assign the labelled array as: NW -> 1, W -> 2, Solid -> 0
        #seg_im_ero_label_show = seg_im_ero_label.copy()
        #seg_im_ero_label_show[seg_im_ero_label_show !=1] = 2
        #seg_im_ero_label_show[np.logical_not(seg_image_2d)]=0
        # ******************** End: for check **************************#
        
        seg_im_ero_label = np.zeros_like(seg_im_ero_label_temp,dtype=bool)
        for labels in label_check:
            seg_im_ero_label = np.logical_or(seg_im_ero_label,seg_im_ero_label_temp==labels)
        seg_im_ero_label = seg_im_ero_label.astype(np.uint8)
        
        # Step 3: perform dilation on the labelled pore space 
        seg_im_ero_label_dil = morphology.binary_dilation(seg_im_ero_label,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of 'seg_im_ero_label_dil' is 'bool'
        seg_im_ero_label_dil = seg_im_ero_label_dil.astype(np.uint8)
        seg_im_ero_label_dil[np.logical_not(seg_im_ero_label_dil)]=2
        seg_im_ero_label_dil[np.logical_not(seg_image)]=0

        Sw = (seg_im_ero_label_dil==2).sum()/pore_vol
    #end if
    return (seg_im_ero_label_dil,Sw)
#end def    

def generate_morph_drain_curv(seg_image_input,R_critical):
    # The method for this function follows Hilper & Miller AWR(2001)
    # 1. Perform erosion for the pore space with radius of R_critical
    # 2. Label the eroded pore space, and leave only the pore space that is still 
    #    connected with the non-wetting phase reservoir
    # 3. Perform the dilation for the labelled pore space with radius of R_critical
    # ****************************************************************************
    # Currently I am provided with a 3D SignDist image which has positive values
    # in the pore space and 0.0 in the solid phase.
    # ****************************************************************************
    if seg_image_input.ndim == 2:

        pore_vol = 1.0*(seg_image_input>0.0).sum()
        radius = R_critical
        print 'Morphological Drainage: processing critical radius: '+str(radius)+' now......'

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)

        # Step 1.2: Perform erosion on the pore space
        # NOTE: the dtype of 'seg_im_ero' is 'bool'
        seg_im_ero = morphology.binary_erosion(seg_image_input>0.0,structure=circle,border_value=1)
        # NOTE: 'border_value' for erosion should be 'True'

        # Step 2: Label the eroded pore space
        # NOTE: Assume the NW phase reservoir is at the first layer of the domain
        #       i.e. at seg_image[0,:] - adjust it if this does not suit your need
        # For erosion, assume that diagonals are not considered
        seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(2,1))
        #seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(2,2))
        # NOTE: Here I assume the inlet is at the first layer of the array's axis=2 (i.e. domain[0,:,:])\
        #       You can always change to any other layers as the inlet for this drainage.
        label_check = seg_im_ero_label_temp[0,seg_im_ero_label_temp[0,:]!=0]
        label_check = np.unique(label_check)

        # NOTE the following lines are only for your to check things
        # ******************** For check *******************************#
        # It assign the labelled array as: NW -> 1, W -> 2, Solid -> 0
        #seg_im_ero_label_show = seg_im_ero_label.copy()
        #seg_im_ero_label_show[seg_im_ero_label_show !=1] = 2
        #seg_im_ero_label_show[np.logical_not(seg_image_2d)]=0
        # ******************** End: for check **************************#
        
        seg_im_ero_label = np.zeros_like(seg_im_ero_label_temp,dtype=bool)
        for labels in label_check:
            seg_im_ero_label = np.logical_or(seg_im_ero_label,seg_im_ero_label_temp==labels)
        #seg_im_ero_label = seg_im_ero_label.astype(np.uint8)
        
        # Step 3: perform dilation on the labelled pore space 
        seg_im_ero_label_dil = morphology.binary_dilation(seg_im_ero_label,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of 'seg_im_ero_label_dil' is 'bool'
        seg_im_ero_label_dil[seg_image_input<=0.0]=False
        
        Sw = 1.0 - seg_im_ero_label_dil.sum()/pore_vol
    else:

        pore_vol = 1.0*(seg_image_input>0.0).sum()
        radius = R_critical
        print 'Morphological Drainage: processing critical radius: '+str(radius)+' now......'

        # Step 1.1: Create structuring element
        domain_size = int(np.rint(radius*2)+2)
        grid = np.indices((domain_size,domain_size,domain_size))
        mk_circle = (grid[0]-domain_size/2)**2 + (grid[1]-domain_size/2)**2 + (grid[2]-domain_size/2)**2 <= radius**2
        circle = np.zeros((domain_size,domain_size,domain_size),dtype=np.uint8)
        circle[mk_circle]=1
        circle = extract_shape(circle).astype(bool)

        # Step 1.2: Perform erosion on the pore space
        # NOTE: the dtype of 'seg_im_ero' is 'bool'
        seg_im_ero = morphology.binary_erosion(seg_image_input>0.0,structure=circle,border_value=1)
        # NOTE: 'border_value' for erosion should be 'True'

        # Step 2: Label the eroded pore space
        # NOTE: Assume the NW phase reservoir is at the first layer of the domain
        #       i.e. at seg_image[0,:] - adjust it if this does not suit your need
        # For erosion, assume that diagonals are not considered
        seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(3,1))
        #seg_im_ero_label_temp,num_features = measurements.label(seg_im_ero,structure=morphology.generate_binary_structure(3,3))
        # NOTE: Here I assume the inlet is at the first layer of the array's axis=2 (i.e. domain[0,:,:])\
        #       You can always change to any other layers as the inlet for this drainage.
        label_check = seg_im_ero_label_temp[0,seg_im_ero_label_temp[0,:]!=0]
        label_check = np.unique(label_check)

        # NOTE the following lines are only for your to check things
        # ******************** For check *******************************#
        # It assign the labelled array as: NW -> 1, W -> 2, Solid -> 0
        #seg_im_ero_label_show = seg_im_ero_label.copy()
        #seg_im_ero_label_show[seg_im_ero_label_show !=1] = 2
        #seg_im_ero_label_show[np.logical_not(seg_image_2d)]=0
        # ******************** End: for check **************************#
        
        seg_im_ero_label = np.zeros_like(seg_im_ero_label_temp,dtype=bool)
        for labels in label_check:
            seg_im_ero_label = np.logical_or(seg_im_ero_label,seg_im_ero_label_temp==labels)
        #seg_im_ero_label = seg_im_ero_label.astype(np.uint8)
        
        # Step 3: perform dilation on the labelled pore space 
        seg_im_ero_label_dil = morphology.binary_dilation(seg_im_ero_label,structure=circle,border_value=0)
        # NOTE: 'border_value' for dilation should be 'False'
        # NOTE: the dtype of 'seg_im_ero_label_dil' is 'bool'
        seg_im_ero_label_dil[seg_image_input<=0.0]=False
        
        Sw = 1.0 - seg_im_ero_label_dil.sum()/pore_vol
    #end if 
    return Sw
#end def    

def extract_shape(domain):
    if domain.ndim == 3:
        where_tube = np.where(domain)
        zStart = where_tube[0].min()
        zEnd   = where_tube[0].max()
        yStart = where_tube[1].min()
        yEnd   = where_tube[1].max()
        xStart = where_tube[2].min()
        xEnd   = where_tube[2].max()
        domain_seg = domain[zStart:zEnd+1,yStart:yEnd+1,xStart:xEnd+1].copy()
        # IMPORTANT: if you have "domain_seg = domain[yStart:yEnd+1,xStart:xEnd+1]"
        #            then "domain_seg" is only a view of domain, and later on you have
        #            any changes on your "domain_seg", the "domain" will also be changed
        #            correspondingly, which might introduce unexpected conflicts and errors
    else: # domain.ndim == 2
        where_tube = np.where(domain)
        yStart = where_tube[0].min()
        yEnd   = where_tube[0].max()
        xStart = where_tube[1].min()
        xEnd   = where_tube[1].max()
        domain_seg = domain[yStart:yEnd+1,xStart:xEnd+1].copy()
    return domain_seg
#end def





