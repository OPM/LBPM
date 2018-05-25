#!/usr/bin/env python

import numpy as np
import sys

# Check if there is a proper command line argument
if len(sys.argv) !=3:
    sys.stderr.write('Usage: ' + sys.argv[0] + ' pores.csv Sw_list\n')
    sys.exit()
# end if

# Read 'pores.csv' and a list of Sw
pores = np.genfromtxt(sys.argv[1])
Sw    = np.genfromtxt(sys.argv[2])
#NOTE: 'pores.csv' has a layout of : [x,y,z,pore_radius]
pores = pores[:,-1]

# Calculate the percentile according to the Sw list
if Sw.max()<=1.0 and Sw.min()>=0.0:
    radii_init_imbib = np.percentile(pores,100.0*Sw.ravel())
else:
    print 'Error: the list of Sw should have values 0.0 - 1.0 !'
    sys.exit()
#end if
radii_init_imbib.shape = (radii_init_imbib.size,1)

# Write out initial guess for imbibition
output_name = 'radius.csv'
print '------ Initial radii for morphological imbibition is written to the file: '+output_name+' ------'
np.savetxt(output_name,radii_init_imbib)


