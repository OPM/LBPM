#!/usr/bin/env python

import sys
import numpy as np
from dist_func_utils import * 
from glob import glob

# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('Usage: ' + sys.argv[0] + ' <Domain.in>\n')
    sys.exit()
# end if

# Read 'Domain.in' to obtain the size of 'SignDist.xxxxx'
f = open(sys.argv[1],'r')
lines = f.readlines()
nx,ny,nz = np.fromstring(lines[1].splitlines()[0],dtype = np.int32,sep=' ')
f.close()
#nx = ny = nz = 128
nx+=2 # consider the halo layer
ny+=2
nz+=2

base_name = 'SignDist.'
file_input_group = glob(base_name+'*')

# Prepare output file: 'pores_xxxxx.csv'
output_data_name = 'pores_'
output_data_format = '.csv'

# Process all imported experimental images
if not file_input_group:
    print 'Error: Input files cannot be found ! '
else:
    for ii in range(len(file_input_group)):
        file_input_single = file_input_group[ii]
        file_input_single_idx = file_input_single[file_input_single.find(base_name)+len(base_name):]
        print '--- Get pore size information for '+file_input_single+' now ---'
        dist = np.fromfile(file_input_single,dtype = np.float64)
        dist_on_skel,connect_stats = detect_intersection_point_with_cluster_filter(dist,nx,ny,nz)
        [z_skel,y_skel,x_skel] = np.where(dist_on_skel>0.0)
        output_data = np.column_stack((x_skel,y_skel,z_skel,dist_on_skel[dist_on_skel>0.0]))
        print '--- Save pore size csv file ---'
        np.savetxt(output_data_name+file_input_single_idx+output_data_format,output_data,delimiter=' ')

    #end for
#end if    









