
import numpy as np
from glob import glob
from lbpm_solid_coordinate_number_utils import *
import matplotlib.pyplot as plt

#NOTE: I could have read the 'Domain.in' files to read the information about the subdomain
#      but let's assume for now all the subdomains are strictly cubic

# set the name of the full domain in *nc format
domain_file_name = "domain1_256.nc"
stats_OUT = solid_coord_fulldomain(domain_file_name)
Aws = cal_Aws_fulldomain(domain_file_name)*1.0 # Convert Aws to a float number
stats_OUT = np.vstack((stats_OUT,np.array([99,Aws]))) #The code 99 is dummy - I just want to attach the Aws to the stats_OUT data
np.savetxt(domain_file_name[:-len(".nc")]+"_stats.txt",stats_OUT)

# Trial plot
plt.figure(1)
plt.semilogy(stats_OUT[1:-1,0],stats_OUT[1:-1,1]/stats_OUT[-1,-1],'ro-')
plt.ylabel('Partial Aws / Total Aws')
plt.xlabel('Number of solid neighbours')
plt.grid(True)
plt.show()



#TODO: make the routine that can analyse individual subdomains and agglomerate the itemfreq data from all subdomains
## Load the ID field "ID.00*"
#ID_prefix="ID."
#halo_layer = 1 # the length of the halo layer
#id_group = glob(ID_prefix+'*')
#id_group.sort()
#
#
#if not id_group:
#    print 'Error: No data files: '+id_group
#else:
#    for ii in range(len(id_group)):
#        print '**Info: Read data file: '+id_group[ii]
#        print "**Info: Start analysing the solid coordinate number......"
#        # call function here
        


