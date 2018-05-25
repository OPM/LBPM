
import numpy as np
from scipy import stats
import scipy.ndimage.morphology as morphology
from netCDF4 import Dataset

def read_NetCDF_file_py(file_name,data_name):

    ncfile = Dataset(file_name,mode='r')
    dataOUT = ncfile.variables[data_name][:]
    
    ncfile.close()
    return dataOUT
#end def

def solid_coord_fulldomain(id_field_file):
    # 'id_field_file' is the file name of the full ID field in *nc format
    # NOTE: This input *nc file is not the raw CT file (which has 1->solid; 0->fluid)
    # Assume the fulldomain is raw, i.e. no postprocessed reservoir layers or porous plate
    # Assume the node type convention for LBPM-WIA simulations, which says
    #        id = 0 -> solid phase
    #        id = 1 -> non-wetting phase
    #        id = 2 -> wetting phase
    # -------------------------------------
    print "**Info: Load the image file: "+id_field_file
    domain = read_NetCDF_file_py(id_field_file,'segmented')
    print "**Info: Start analysing the solid coordinate number ......"
    domain = np.logical_not(domain) # Now 1 -> solid nodes; 0 -> fluid nodes
    domain = domain.astype(np.int8)
    
    # Define the D3Q19 lattice
    cx=np.array([0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0, 0])
    cy=np.array([0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1,-1])
    cz=np.array([0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1, 1])
    
    z_axis = 2
    y_axis = 1
    x_axis = 0
    domain_temp = np.zeros_like(domain)
    for idx in np.arange(1,19):
        domain_temp += np.roll(np.roll(np.roll(domain,cx[idx],axis=x_axis),cy[idx],axis=y_axis),cz[idx],axis=z_axis)
    #end for
    domain_temp = domain_temp[domain==0] # only extract the coordinate number for pore space nodes 
                                         # NOTE that we have 0 -> fluid nodes and 1 -> solid nodes in domain
    return stats.itemfreq(domain_temp)
#end def

def cal_Aws_fulldomain(id_field_file):
    # Calculate the (total) fluid-solid interfacial area for an input porous medium 'Aws'
    # 'id_field_file' is the file name of the full ID field in *nc format
    # NOTE: This input *nc file is not the raw CT file (which has 1->solid; 0->fluid)
    # Assume the fulldomain is raw, i.e. no postprocessed reservoir layers or porous plate
    # Assume the node type convention for LBPM-WIA simulations, which says
    #        id = 0 -> solid phase
    #        id = 1 -> non-wetting phase
    #        id = 2 -> wetting phase
    # -------------------------------------
    print "**Info: Load the image file: "+id_field_file
    domain = read_NetCDF_file_py(id_field_file,'segmented')
    domain[domain>0]=1 # in case there are some nodes with id=2
    print "**Info: Start calculating the fluid-solid interfacial area ......"
    
    # Generate a D3Q19 structure unit cell
    D3Q19=morphology.generate_binary_structure(3,2)
    domain_dil=morphology.binary_dilation(domain,structure=D3Q19,border_value=0).astype(np.int8)
    #NOTE: It is important to have 'border_value' set as 0 for morphological dilation !
    domain_dil = domain_dil - domain
    if (domain_dil<0).sum() > 0:
        print "**Error: The domain for calculating the fluid-solid interfacial area is not set properly !"
    #end if
    Aws=domain_dil.sum()
    return Aws
#end def

#TODO: implement a similar function to analyse only the subdomain 
#      and find a way to concatenate the itemfeq data from all subdomains
def solid_coord_subdomain(id_field_file):
    #NOTE: This function is not done yet *******************************************
    # 'id_field_file' is the file name of the ID field, i.e. "ID.00*"
    # Assume the input doamin is a subdomain segmented by lbpm_segment_decomp.cpp
    #        thus there is an extra halo layer wrapping the subdomain
    # Assume the subdomain is strictly cubic
    # Assume the subdomain is raw, i.e. originally segmented CT image with no reservoir layers or porous plate
    # Assume the node type convention for LBPM-WIA simulations, which says
    #        id = 0 -> solid phase
    #        id = 1 -> non-wetting phase
    #        id = 3 -> wetting phase
    # -------------------------------------
    domain = np.fromfile(id_field_file,dtype=np.int8)
    l_domain = int(pow(domain.size,1.0/3.0))
    domain.shape=(l_domain,l_domain,l_domain)
    domain = np.logical_not(domain) # Now 1 -> solid nodes; 0 -> fluid nodes
    domain = domain.astype(np.int8)

    # Define the D3Q19 lattice
    cx=np.array([0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0, 0])
    cy=np.array([0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1,-1])
    cz=np.array([0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1, 1])
    
    z_axis = 2
    y_axis = 1
    x_axis = 0
    domain_temp = np.zeros_like(domain)
    for idx in np.arange(1,19):
        domain_temp += np.roll(np.roll(np.roll(domain,cx[idx],axis=x_axis),cy[idx],axis=y_axis),cz[idx],axis=z_axis)
    #end for

#end def




