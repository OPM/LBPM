#!/usr/bin/env python

import numpy as np
from scipy import ndimage
from scipy import spatial
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.stats as stats
from skimage.morphology import medial_axis,skeletonize_3d

def detect_intersection_point_with_cluster_filter(input_data,lx,ly,lz=0):
    # Input is a signed distance data file such as 'SignDist.xxxxx' 
    if lz > 0: # i.e. a 3D input
        # Calculate the medial axis of the segmented image
        dist = input_data.copy() 
        dist.shape = (lz,ly,lx)
        skel = skeletonize_3d(dist>0.0)
        dist_on_skel = skel*dist
        [z_skel,y_skel,x_skel]=np.where(skel)
        
        # Building two search trees is potentially good for 3D image
        grid = np.indices(dist.shape)
        grid_points = zip(grid[0].ravel(),grid[1].ravel(),grid[2].ravel())
        tree_grid = spatial.cKDTree(grid_points)
        points_for_search = zip(z_skel,y_skel,x_skel)
        tree_points_for_search = spatial.cKDTree(points_for_search)
        neighbor_all = tree_points_for_search.query_ball_tree(tree_grid,np.sqrt(3.0))

        idx_glb_table = np.ravel_multi_index([z_skel,y_skel,x_skel],skel.shape)
        idx_glb_candidate = np.empty((0,),dtype=int)
        for k,idx_glb_neighbor in enumerate(neighbor_all):
            if k%4000==0:
                print 'Search for intersection points: '+str(k)+'/'+str(len(neighbor_all)-1)
            #end if
            mask = np.in1d(idx_glb_neighbor,idx_glb_table,assume_unique=True)
            idx_glb_candidate = np.hstack((idx_glb_candidate,mask.sum()))
        #end for   

        # Statistics: number of neighbors for each voxel on the medial axis
        connectivity_stats = stats.itemfreq(idx_glb_candidate)
        # 'connectivity_stats' has the following format:
        # number_of_neighbors_for_each_voxel :: corresponding_number_of_voxels
        # array([[  1,  41],
        #        [  2, 143],
        #        [  3, 185],
        #        [  4,   5]])
        #
        # If a voxel is indentified as an intersection point, the number of its
        # neighboring voxels should be greater than 'benchmark'
        benchmark = connectivity_stats[np.argmax(connectivity_stats[:,1]),0]
        # Update the medial axis and 'dist_on_skel'
        skel[:] = False
        skel[np.unravel_index(idx_glb_table[idx_glb_candidate > benchmark],skel.shape)] = True
    else: # i.e. 2D input
        # Calculate the medial axis of the segmented image
        dist = input_data.copy() 
        dist.shape = (ly,lx)
        skel = medial_axis(dist>0.0)
        dist_on_skel = skel*dist
        [y_skel,x_skel]=np.where(skel)
        
        # Building two search trees is potentially good for 3D image
        grid = np.indices(dist.shape)
        grid_points = zip(grid[0].ravel(),grid[1].ravel())
        tree_grid = spatial.cKDTree(grid_points)
        points_for_search = zip(y_skel,x_skel)
        tree_points_for_search = spatial.cKDTree(points_for_search)
        neighbor_all = tree_points_for_search.query_ball_tree(tree_grid,np.sqrt(2.0))

        idx_glb_table = np.ravel_multi_index([y_skel,x_skel],skel.shape)
        idx_glb_candidate = np.empty((0,),dtype=int)
        for k,idx_glb_neighbor in enumerate(neighbor_all):
            if k%4000==0:
                print 'Search for intersection points: '+str(k)+'/'+str(len(neighbor_all)-1)
            #end if
            mask = np.in1d(idx_glb_neighbor,idx_glb_table,assume_unique=True)
            idx_glb_candidate = np.hstack((idx_glb_candidate,mask.sum()))
        #end for   

        # Statistics: number of neighbors for each voxel on the medial axis
        connectivity_stats = stats.itemfreq(idx_glb_candidate)
        # 'connectivity_stats' has the following format:
        # number_of_neighbors_for_each_voxel :: corresponding_number_of_voxels
        # array([[  1,  41],
        #        [  2, 143],
        #        [  3, 185],
        #        [  4,   5]])
        #
        # If a voxel is indentified as an intersection point, the number of its
        # neighboring voxels should be greater than 'benchmark'
        benchmark = connectivity_stats[np.argmax(connectivity_stats[:,1]),0]
        # Update the medial axis and 'dist_on_skel'
        skel[:] = False
        skel[np.unravel_index(idx_glb_table[idx_glb_candidate > benchmark],skel.shape)] = True
    #end if
    return (_Filter_cluster_close_points(skel*dist),connectivity_stats)
#end def    

def _Filter_cluster_close_points(arr):
    # 'arr' is a 2D/3D signed distance data file
    # Return an 'arr' with nearest neighboring points clustered
    if arr.ndim == 2:
        [y_idx,x_idx] = np.where(arr>0.0)
        idx_glb = np.ravel_multi_index([y_idx,x_idx],arr.shape)
        grid = np.indices(arr.shape)
        
        dist = arr.ravel()[idx_glb]
        candidate = np.empty((0,),dtype=int)
        # TODO: use search tree to do this !
        for k in np.arange(idx_glb.size):
            if k%200==0:
                print 'Clustering close points: '+str(k)+'/'+str(idx_glb.size-1)
            #end if
            mask_circle = (grid[0]-y_idx[k])**2 + (grid[1]-x_idx[k])**2<=dist[k]**2
            idx_glb_circle =np.ravel_multi_index(np.where(mask_circle),mask_circle.shape)
            mask = np.in1d(idx_glb,idx_glb_circle,assume_unique=True)
            candidate = np.hstack((candidate,idx_glb[mask][np.argmax(dist[mask])])) 
        #end for    
        candidate = np.unique(candidate)
        mask = np.in1d(idx_glb,candidate,assume_unique=True)
        output_glb_idx = idx_glb[mask]
        output = np.zeros_like(arr.ravel())
        output[output_glb_idx] = arr.ravel()[output_glb_idx]
        output.shape = arr.shape
    else:
        [z_idx,y_idx,x_idx] = np.where(arr>0.0)
        idx_glb = np.ravel_multi_index([z_idx,y_idx,x_idx],arr.shape)
        grid = np.indices(arr.shape)
        
        dist = arr.ravel()[idx_glb]
        candidate = np.empty((0,),dtype=int)
        # TODO: use search tree to do this !
        for k in np.arange(idx_glb.size):
            if k%200==0:
                print 'Clustering close points: '+str(k)+'/'+str(idx_glb.size-1)
            #end if
            mask_circle = (grid[0]-z_idx[k])**2 + (grid[1]-y_idx[k])**2 + (grid[2]-x_idx[k])**2<=dist[k]**2
            idx_glb_circle =np.ravel_multi_index(np.where(mask_circle),mask_circle.shape)
            mask = np.in1d(idx_glb,idx_glb_circle,assume_unique=True)
            candidate = np.hstack((candidate,idx_glb[mask][np.argmax(dist[mask])])) 
        #end for    
        candidate = np.unique(candidate)
        mask = np.in1d(idx_glb,candidate,assume_unique=True)
        output_glb_idx = idx_glb[mask]
        output = np.zeros_like(arr.ravel())
        output[output_glb_idx] = arr.ravel()[output_glb_idx]
        output.shape = arr.shape
    #end if
    return output
#end def    

def get_Dist(f):
    return  ndimage.distance_transform_edt(f)
#end def    

def detect_local_maxima(arr,medial_axis_arr,patch_size=3):
    local_max = _find_local_maxima(arr,patch_size)
    background_value = 1e5 # for denoising in finding local minima
    arr2 = arr.copy()
    arr2[np.logical_not(medial_axis_arr)]=background_value
    ocal_min = _find_local_minima(arr2,background_val=background_value)
    local_min = np.logical_and(medial_axis_arr,local_min) #Correct min_indices with Medial axis
    local_max_reduced = np.logical_xor(local_max,local_min)
    local_max = np.logical_and(local_max,local_max_reduced)
    local_max = np.logical_and(medial_axis_arr,local_max)#Correct max_indices with Medial axis
    return local_max
#end def    

def detect_local_maxima_with_cluster_filter(arr,medial_axis_arr,patch_size=3):
    # apply the cluster filetering only once
    local_max = detect_local_maxima(arr,medial_axis_arr,patch_size=patch_size)
    arr2 = arr.copy()
    arr2[np.logical_not(local_max)]=0.0
    return _Filter_cluster_close_points(arr2)
#end def    

def detect_local_maxima_with_cluster_filter_loop(arr,medial_axis_arr,patch_start=3,patch_stop=9):
    # Implement multiple search for local maxima (with cluster filtering)
    arr1 = detect_local_maxima_with_cluster_filter(arr,medial_axis_arr,patch_size=patch_start) 
    [y_idx_1,x_idx_1] = np.where(arr1>0.0)
    arr2 = detect_local_maxima_with_cluster_filter(arr1,medial_axis_arr,patch_size=patch_start+2) 
    [y_idx_2,x_idx_2] = np.where(arr2>0.0)
    if (y_idx_1.size>y_idx_2.size):
        counter = 2 # Record how many times is the filtering applied
        patch = np.arange((patch_start+2)+2,patch_stop+2,2)
        for k in patch:
            
            y_idx_1 = y_idx_2.copy()
            x_idx_1 = x_idx_2.copy()

            arr1 = detect_local_maxima_with_cluster_filter(arr2,medial_axis_arr,patch_size=k) 
            [y_idx_2,x_idx_2] = np.where(arr1>0.0)
            
            arr2 = arr1.copy()
            
            counter = counter + 1
            if not (y_idx_1.size > y_idx_2.size):
                print 'Note: Local maxima are already found before the patch_stop is reached.'
                print '      All local maxima are found at the patch size of '+str(k)+' !'
                break;
            #end if
        #end for
    else:
        print 'Note: All local maxima are found at the patch size of '+str(patch_start)+' !'
        return arr2
    #end if
    return arr2
#end def  

#def Filter_cluster(arr):
#    return _Filter_cluster_close_points(arr) 
##end def
def _find_local_maxima(arr,patch_size):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    #neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    #neighborhood = np.ones((patch_size,patch_size),dtype=bool)
    # apply the local maximum filter; all locations of maximum value 
    # in their neighborhood are set to 1
    local_max = (filters.maximum_filter(arr, size=patch_size,mode='constant',cval=0)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    if arr.ndim == 3:
        neighborhood = np.ones((patch_size,patch_size,patch_size),dtype=bool)
    elif arr.ndim==2:
        neighborhood = np.ones((patch_size,patch_size),dtype=bool)
    else:
        neighborhood = np.ones((patch_size,),dtype=bool)
    #end if
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_maxima = local_max - eroded_background
    return detected_maxima

def _find_local_minima(arr,patch_size=3,background_val=1e5):
    local_min = (filters.minimum_filter(arr, size=patch_size,mode='constant',cval=background_val)==arr)
    background = (arr==background_val)
    if arr.ndim == 3:
        neighborhood = np.ones((patch_size,patch_size,patch_size),dtype=bool)
    elif arr.ndim==2:
        neighborhood = np.ones((patch_size,patch_size),dtype=bool)
    else:
        neighborhood = np.ones((patch_size,),dtype=bool)
    #end if
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    detected_minima = local_min - eroded_background
    return detected_minima
#end def



