#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from glob import glob
import sys
#from dist_func_utils import *
import scipy.stats as stats
import scipy.ndimage.morphology as morphology
import scipy.ndimage.measurements as measurements
from morphological_analysis_utils import *

# Check if there is a proper command line argument
#if len(sys.argv) !=2:
#    sys.stderr.write('Usage: ' + sys.argv[0] + ' <Base name of the segmented image data>\n')
#    sys.exit()
# end if

base_name = 'MicroModel' # Give the name of the input parameter file 
#base_name = 'SignDist' # Give the name of the input parameter file 
#base_name = sys.argv[1] # Give the name of the input parameter file 
SegImage_format = '.raw'

# Load a group of segmented image data files with 'base_name' in the names of files 
# e.g. 'Micromodel_1_segmented.raw' etc.
file_input_group = glob('*'+base_name+'*'+SegImage_format) # need to match the image data format

# Dimensions for segmented image
lz = 500
ly = 500
lx = 12
#lz = 200
#ly = 200
#lx = 200

R_critical = 4.5

# Process all imported experimental images
if not file_input_group:
    print 'Error: Input files cannot be found ! '
else:
    for ii in range(len(file_input_group)):
        file_input_single = file_input_group[ii]
        print "Analyzing segmented image data "+file_input_single+" now..."

        # The imported data is a double-precision signed distance map
        image_raw = np.fromfile(file_input_single,dtype=np.uint8)
        #image_raw = np.fromfile(file_input_single,dtype=np.float64)
        image_raw.shape=(lz,ly,lx)
        
        edt_dist = morphology.distance_transform_edt(image_raw)
        #edt_dist = edt_dist[:,:,5]
        # wrap the medium with one layer of solid phase 
        #image_raw = np.lib.pad(image_raw,((0,0),(1,1),(1,1)),'constant',constant_values=0.0)

        drain_config,Sw_drain = generate_morph_drain_config(edt_dist,R_critical)
        imbib_config,Sw_imbib = generate_morph_imbib_config(edt_dist,R_critical)
        #drain_config,Sw_drain = generate_morph_drain_config(image_raw,R_critical)
        #imbib_config,Sw_imbib = generate_morph_imbib_config(image_raw,R_critical)
        #Sw_drain = generate_morph_drain_curv_3D(image_raw,R_critical)
        #imbib_config,Sw_imbib = generate_morph_imbib_config2_3D(image_raw,R_critical)
        #Sw_imbib = generate_morph_imbib_curv2_3D(image_raw,R_critical)
        print 'Rcrit: '+str(R_critical)+',  Drain_Sw = '+str(Sw_drain)
        print 'Rcrit: '+str(R_critical)+',  Imbib_Sw = '+str(Sw_imbib)

        # Save the configuration files
#        drain_config.tofile('Drain_config_Rcrit_'+str(R_critical)+'.raw')
#        imbib_config.tofile('Imbib_config_Rcrit_'+str(R_critical)+'.raw')

        # Load the saved data
#        drain_config = np.fromfile('Drain_config_Rcrit_'+str(R_critical)+'.raw',dtype=np.uint8)
#        imbib_config = np.fromfile('Imbib_config_Rcrit_'+str(R_critical)+'.raw',dtype=np.uint8)
#        drain_config.shape = (lz,ly,lx)
#        imbib_config.shape = (lz,ly,lx)


        plt.figure(1)
        plt.subplot(1,2,1)
        plt.title('Drainage: Rcrit='+str(R_critical)+' Sw='+str(Sw_drain))
        plt.pcolormesh(drain_config[:,:,lx/2],cmap='hot')
        #plt.pcolormesh(drain_config[20,:,:],cmap='hot')
        plt.axis('equal')
        plt.grid(True)

        plt.subplot(1,2,2)
        plt.title('Imbibition: Rcrit='+str(R_critical)+' Sw='+str(Sw_imbib))
        plt.pcolormesh(imbib_config[:,:,lx/2],cmap='hot')
        #plt.pcolormesh(imbib_config[20,:,:],cmap='hot')
        plt.axis('equal')
        plt.grid(True)

#        plt.figure(1)
#        plt.subplot(1,2,1)
#        plt.title('Drainage: Rcrit='+str(R_critical)+' Sw='+str(Sw_drain))
#        plt.pcolormesh(drain_config,cmap='hot')
#        plt.axis('equal')
#        plt.grid(True)
#
#        plt.subplot(1,2,2)
#        plt.title('Imbibition: Rcrit='+str(R_critical)+' Sw='+str(Sw_imbib))
#        plt.pcolormesh(imbib_config,cmap='hot')
#        plt.axis('equal')
#        plt.grid(True)

        plt.figure(2)
        plt.plot(Sw_drain,1.0/R_critical,'ro',markersize=6,label='Drainage')
        plt.plot(Sw_imbib,1.0/R_critical,'b^',markersize=6,label='Imbibition')
        plt.legend(loc='best')
        plt.grid(True)
        plt.xlim([0,1.0])
        plt.ylim([0,2.0])

        

        plt.show()
    #end if
#end for    


