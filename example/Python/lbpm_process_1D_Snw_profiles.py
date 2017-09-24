
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from glob import glob

def read_NetCDF_file_py(file_name,data_name):

    ncfile = Dataset(file_name,mode='r')
    dataOUT = ncfile.variables[data_name][:]
    
    ncfile.close()
    return dataOUT
#end def

# Load the domain
domain = read_NetCDF_file_py('domain1_256_reduced_lbpm.nc','segmented')
vol_pore_1D = np.sum(np.sum(domain,axis=2),axis=1)*1.0
# Set the files & names of the phase data 
data_group_prefix="OutputData_vis*.nc"
data_group = glob(data_group_prefix)
#data_group.sort()
data_name="Phase"

output_data = np.empty((0,),dtype=np.float64)
if not data_group:
    print "**Error: No input file group: "+data_group_prefix
else:
    for ii in np.arange(len(data_group)):
        data_single = data_group[ii]
        print "**Info: Start processing file "+data_single+" now..."
        time_step = int(data_single[data_single.find("_vis")+len("_vis"):data_single.find(".nc")])
        phase = read_NetCDF_file_py(data_group[ii],data_name)
        domain_nw = (phase>0.0).astype(np.int8)
        vol_nw_1D = np.sum(np.sum(domain_nw,axis=2),axis=1)*1.0
        Snw_1D = vol_nw_1D/vol_pore_1D
        Snw_1D[Snw_1D==0.0] = 1.0e-10
        output_data = np.hstack((output_data,np.hstack((time_step,Snw_1D))))
        

    #end for
#end if    
output_data.shape = (len(data_group),output_data.size/len(data_group))
output_data = output_data[output_data[:,0].argsort(),:] #sort output data according to time_step

print "**Info: Save the 1D Snw data now."
np.savetxt('1D_Snw_data.txt',output_data)


