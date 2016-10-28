#!/usr/bin/env python

import numpy as np
import sys
import os
import shutil
from preprocess_utils import *

# 1. this script take one input argument: experiment.csv
# 2. Users need to put two things into this script
#    2.1 surface tension in physical units
#    2.2 physical depth
# 3. Read csv + convert to LBM information
# 4. write Color.in,-> this needs info from LBM pressure
# 5. write Segmented.in -> this needs names of segmented image data
# 6. write Domain.in -> this needs info on how users want to decompose the domain

# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('Usage: ' + sys.argv[0] + ' <Input experiment file>\n')
    sys.exit()
# end if

#experiment_file = 'experiment.csv' # Give the name of the experiment file (e.g. *.csv)
experiment_file = sys.argv[1] # Give the name of the experiment file (e.g. *.csv)
seg_image_data_suffix = '_segmented.raw'# suffix of the segmented image data
                                        # (should be consistent with what's in the segmenting script)
                                        # TODO: It'd be good to find a better way to do this
image_format = '.tiff'
process_name = 'drainage'

#### Users need to put information here ####
ift = 24.0   # dyne/cm
Depth = 8.8  # micron

# A list of all default values in 'Color.in'
#    Para['tau']         = 0.7
#    Para['alpha']       = 0.005
#    Para['beta']        = 0.95
#    Para['phi_solid']   = -1.0
#    Para['saturation']  = 0.0
#    Para['Fx']          = 0.0
#    Para['Fy']          = 0.0
#    Para['Fz']          = 0.0
#    Para['Restart']     = 0
#    Para['pBC']         = 1
#    Para['din']         = 1.001
#    Para['dout']        = 0.999
#    Para['maxtime']     = 100005
#    Para['interval']    = 2000
#    Para['tolerance']   = 1e-5

# ***** Update any variables in 'Color.in', using names given in the key ***** # 
alpha = 0.01

# **************************************************************************** #

# A list of all default values in 'Domain.in'
#    Para['nprocx']       = 1
#    Para['nprocy']       = 2
#    Para['nprocz']       = 2
#    Para['nx']           = 1
#    Para['ny']           = 2
#    Para['nz']           = 2
#    Para['nspheres']     = 0 # deprecated
#    Para['Lx']           = 10
#    Para['Ly']           = 500
#    Para['Lz']           = 500

# ***** Update any variables in 'Domain.in', using names given in the key ***** # 


# ***************************************************************************** #

# A list of all default values in 'Segmented.in'
#    Para['file_name']     = 'Micromodel_1_segmented.raw'
#    Para['Nx']            = 10
#    Para['Ny']            = 500
#    Para['Nz']            = 500
#    Para['xStart']        = 0
#    Para['yStart']        = 0
#    Para['zStart']        = 0

# ***** Update any variables in 'Segmented.in', using names given in the key ***** # 


# ******************************************************************************** #

# Extract key parameters for LBM simulation from the experimental input *.csv file
(Seg_data_name,din,dout)=get_LBM_parameters(experiment_file,seg_image_data_suffix,image_format,\
                                            ift=ift,Depth=Depth)
# Now 'name_for_Segmented_in' should match the name of segmented data files that are already generated

# Write out 'Color.in', 'Domain.in' and 'Segmented.in' files
cwd = os.getcwd()
for k in range(Seg_data_name.size):
    tag = k+1 # tag number for different folders to be created
    dir_name = process_name+'_'+str(tag)
    print "Creating folder : "+dir_name
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    #end if
    
    # Either move the corresponding '*.raw' data file to the folder 'dir_name'
    #os.rename('./'+Seg_data_name[k],'./'+dir_name+'/'+Seg_data_name[k])
    # Or copy the corresponding '*.raw' data file to the folder 'dir_name'
    shutil.copy('./'+Seg_data_name[k],'./'+dir_name)
    
    # Change to the corresponding folder and write all input files 
    os.chdir(dir_name)
    write_Color_in_file(din=din[k],dout=dout[k],alpha=alpha)
    write_Segment_in_file(Seg_data_name=Seg_data_name[k])
    write_Domain_in_file()

    os.chdir(cwd)
#end for    







