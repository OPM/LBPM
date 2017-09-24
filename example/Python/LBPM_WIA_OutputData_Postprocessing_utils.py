#!/usr/bin/env python

#import os
#import sys
import numpy as np
#import matplotlib.pyplot as plt
from glob import glob
from netCDF4 import Dataset

def Writer_NetCDF_multiVar_LBPMWIA(file_name,data_IN,LX,LY,LZ):
    #TODO: add some extra **kwargs so that users can pick a paricular variable name 
    #      and write out what data they want to write out
    #NOTE (LX,LY,LZ) is the size of a single variable
    #NOTE data_IN contains 1D data with the format (Phase, Pressure, BlobID)
    # open a new netCDF file for writing.
    ncfile = Dataset(file_name,'w') 
    # create the output data.
    
    # create the x, y and z dimensions.
    ncfile.createDimension('x',LX)
    ncfile.createDimension('y',LY)
    ncfile.createDimension('z',LZ)
    # create the variable (4 byte integer in this case)
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    Phase    = ncfile.createVariable('Phase',np.float32,('z','y','x'))
    Pressure = ncfile.createVariable('Pressure',np.float32,('z','y','x'))
    BlobID   = ncfile.createVariable('BlobID',np.float32,('z','y','x'))
    Phase[:]    = data_IN[:LX*LY*LZ].reshape((LZ,LY,LX))
    Pressure[:] = data_IN[LX*LY*LZ:2*LX*LY*LZ].reshape((LZ,LY,LX))
    BlobID[:]   = data_IN[2*LX*LY*LZ:].reshape((LZ,LY,LX))
    #data.voxel_unit = 'micrometer'
    #data.voxel_size=5.7
    ncfile.close()
    print '**Info: The *.nc file is written successfully !'
#end def

def Writer_NetCDF_multiVar_Restart_LBPMWIA(file_name,data_IN,LX,LY,LZ):
    #NOTE (LX,LY,LZ) is the size of a single variable
    #NOTE data_IN contains 1D data with the format (Phase, Pressure, BlobID)
    # open a new netCDF file for writing.
    ncfile = Dataset(file_name,'w') 

    # Total grid size
    N=LX*LY*LZ

    # create the x, y and z dimensions.
    ncfile.createDimension('x',LX)
    ncfile.createDimension('y',LY)
    ncfile.createDimension('z',LZ)
    # create the variable (4 byte integer in this case)
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    den_nw=ncfile.createVariable('DensityNW',np.float32,('z','y','x'))
    den_w=ncfile.createVariable('DensityW',np.float32,('z','y','x'))
    f0=ncfile.createVariable('f0',np.float32,('z','y','x'))
    f1=ncfile.createVariable('f1',np.float32,('z','y','x'))
    f2=ncfile.createVariable('f2',np.float32,('z','y','x'))
    f3=ncfile.createVariable('f3',np.float32,('z','y','x'))
    f4=ncfile.createVariable('f4',np.float32,('z','y','x'))
    f5=ncfile.createVariable('f5',np.float32,('z','y','x'))
    f6=ncfile.createVariable('f6',np.float32,('z','y','x'))
    f7=ncfile.createVariable('f7',np.float32,('z','y','x'))
    f8=ncfile.createVariable('f8',np.float32,('z','y','x'))
    f9=ncfile.createVariable('f9',np.float32,('z','y','x'))
    f10=ncfile.createVariable('f10',np.float32,('z','y','x'))
    f11=ncfile.createVariable('f11',np.float32,('z','y','x'))
    f12=ncfile.createVariable('f12',np.float32,('z','y','x'))
    f13=ncfile.createVariable('f13',np.float32,('z','y','x'))
    f14=ncfile.createVariable('f14',np.float32,('z','y','x'))
    f15=ncfile.createVariable('f15',np.float32,('z','y','x'))
    f16=ncfile.createVariable('f16',np.float32,('z','y','x'))
    f17=ncfile.createVariable('f17',np.float32,('z','y','x'))
    f18=ncfile.createVariable('f18',np.float32,('z','y','x'))

    den_nw[:]  = data_IN[   :1*N].reshape(LZ,LY,LX)
    den_w[:]   = data_IN[1*N:2*N].reshape(LZ,LY,LX)
    f0[:]      = data_IN[2*N:3*N].reshape(LZ,LY,LX)
    f1[:]      = data_IN[3*N:4*N].reshape(LZ,LY,LX)
    f2[:]      = data_IN[4*N:5*N].reshape(LZ,LY,LX)
    f3[:]      = data_IN[5*N:6*N].reshape(LZ,LY,LX)
    f4[:]      = data_IN[6*N:7*N].reshape(LZ,LY,LX)
    f5[:]      = data_IN[7*N:8*N].reshape(LZ,LY,LX)
    f6[:]      = data_IN[8*N:9*N].reshape(LZ,LY,LX)
    f7[:]      = data_IN[9*N:10*N].reshape(LZ,LY,LX)
    f8[:]      = data_IN[10*N:11*N].reshape(LZ,LY,LX)
    f9[:]      = data_IN[11*N:12*N].reshape(LZ,LY,LX)
    f10[:]     = data_IN[12*N:13*N].reshape(LZ,LY,LX)
    f11[:]     = data_IN[13*N:14*N].reshape(LZ,LY,LX)
    f12[:]     = data_IN[14*N:15*N].reshape(LZ,LY,LX)
    f13[:]     = data_IN[15*N:16*N].reshape(LZ,LY,LX)
    f14[:]     = data_IN[16*N:17*N].reshape(LZ,LY,LX)
    f15[:]     = data_IN[17*N:18*N].reshape(LZ,LY,LX)
    f16[:]     = data_IN[18*N:19*N].reshape(LZ,LY,LX)
    f17[:]     = data_IN[19*N:20*N].reshape(LZ,LY,LX)
    f18[:]     = data_IN[20*N:21*N].reshape(LZ,LY,LX)
    ncfile.close()
    print '**Info: The *.nc file is written successfully !'
#end def


def Writer_NetCDF_singleVar_LBPMWIA(file_name,data_name,data_IN,LX,LY,LZ):
	
    #NOTE (LX,LY,LZ) is the size of a single variable
    #NOTE data_IN contains 1D data with the format (Phase, Pressure, BlobID)
    # open a new netCDF file for writing.
    ncfile = Dataset(file_name,'w') 
    # create the output data.
	
    # create the x, y and z dimensions.
    ncfile.createDimension('x',LX)
    ncfile.createDimension('y',LY)
    ncfile.createDimension('z',LZ)
    # create the variable (4 byte integer in this case)
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    data    = ncfile.createVariable(data_name,np.float32,('z','y','x'))
    data[:] = data_IN.reshape((LZ,LY,LX))
    #data.voxel_unit = 'micrometer'
    #data.voxel_size=5.7
    ncfile.close()
    print '**Info: The *.nc file is written successfully !'
#end def

def Reconstruct_single_data(data_name_group,proc_config,grid_config):
    #NOTE: This function is not used for reconstructing ID.* and SignDist.*
    #data_name_group is a string, e.g. 'Phase.00*'
    #proc_config is a tuple with (nproc_x,nproc_y,nproc_z), the processor configuration
    #NOTE: grid_config is a tuple with (lx,ly,lz) is the size of a segmented domain
    data_group=glob(data_name_group)
    data_group.sort()
    nproc_x, nproc_y, nproc_z = proc_config
    lx,ly,lz = grid_config
    print '**Info: data reconstruction: '+data_group[0]
    for z_axis in np.arange(nproc_z):
        for y_axis in np.arange(nproc_y):
            for x_axis in np.arange(0,nproc_x,2):
                if nproc_x%2==0: # nproc_x is an even number
                    temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                    temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x2.shape = (lz,ly,lx)
                    #temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]# for ID.* and SignDist.* which have one extra buffer layer around the data
                    #temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                    temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                else: # nproc_x is an odd number
                    if nproc_x==1:
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                        temp_x1.shape = (lz,ly,lx)
                    elif x_axis < nproc_x-2:
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                        temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                        temp_x1.shape = (lz,ly,lx)
                        temp_x2.shape = (lz,ly,lx)
                        #temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                        #temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                        temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                    #end if
                #end if
                if x_axis == 0:
                    ID_x_temp = temp_x1.copy()
                else:
                    ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
                #end if
            #end for
            if nproc_x !=1 and nproc_x%2 != 0:
                temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+nproc_x-1],dtype=np.float64)
                ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
            #end if

            if y_axis == 0:
                ID_y_temp = ID_x_temp.copy()
            else:
                ID_y_temp = np.concatenate((ID_y_temp,ID_x_temp),axis=1)
            #end if
        #end for

        if z_axis == 0:
            ID_z_temp = ID_y_temp.copy()
        else:
            ID_z_temp = np.concatenate((ID_z_temp,ID_y_temp),axis=0)
        #end if
    #end for
    print '**Info: data reconstruction completed for '+data_group[0]
    return ID_z_temp.flatten()
#end def

def Reconstruct_SignDist_data(data_name_group,proc_config,grid_config):
    #NOTE: This function is only used for reconstructing ID.* and SignDist.*
    #data_name_group is a string, e.g. 'ID.00*'
    #proc_config is a tuple with (nproc_x,nproc_y,nproc_z), the processor configuration
    #NOTE: grid_config is a tuple with (lx,ly,lz) is the size of a segmented domain
    data_group=glob(data_name_group)
    data_group.sort()
    nproc_x, nproc_y, nproc_z = proc_config
    lx,ly,lz = grid_config
    print '**Info: data reconstruction: '+data_group[0]
    for z_axis in np.arange(nproc_z):
        for y_axis in np.arange(nproc_y):
            for x_axis in np.arange(0,nproc_x,2):
                if nproc_x%2==0: # nproc_x is an even number
                    temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                    temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x2.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]# for ID.* and SignDist.* which have one extra buffer layer around the data
                    temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                    temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                else: # nproc_x is an odd number
                    if nproc_x==1:
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                        temp_x1.shape = (lz,ly,lx)
                        temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                    elif x_axis < nproc_x-2: # i.e. nproc_x = 3 or 5 or 7 ...
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.float64)
                        temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.float64)
                        temp_x1.shape = (lz,ly,lx)
                        temp_x2.shape = (lz,ly,lx)
                        temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                        temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                        temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                    #end if
                #end if
                if x_axis == 0:
                    ID_x_temp = temp_x1.copy()
                else:
                    ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
                #end if
            #end for
            if nproc_x !=1 and nproc_x%2 != 0:
                temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+nproc_x-1],dtype=np.float64)
                ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
            #end if

            if y_axis == 0:
                ID_y_temp = ID_x_temp.copy()
            else:
                ID_y_temp = np.concatenate((ID_y_temp,ID_x_temp),axis=1)
            #end if
        #end for

        if z_axis == 0:
            ID_z_temp = ID_y_temp.copy()
        else:
            ID_z_temp = np.concatenate((ID_z_temp,ID_y_temp),axis=0)
        #end if
    #end for
    print '**Info: Signed distance data reconstruction completed for '+data_group[0]
    return ID_z_temp.flatten()
#end def

def Reconstruct_ID_data(data_name_group,proc_config,grid_config):
    #NOTE: This function is only used for reconstructing ID.* and SignDist.*
    #data_name_group is a string, e.g. 'ID.00*'
    #proc_config is a tuple with (nproc_x,nproc_y,nproc_z), the processor configuration
    #NOTE: grid_config is a tuple with (lx,ly,lz) is the size of a segmented domain
    data_group=glob(data_name_group)
    data_group.sort()
    nproc_x, nproc_y, nproc_z = proc_config
    lx,ly,lz = grid_config
    print '**Info: data reconstruction: '+data_group[0]
    for z_axis in np.arange(nproc_z):
        for y_axis in np.arange(nproc_y):
            for x_axis in np.arange(0,nproc_x,2):
                if nproc_x%2==0: # nproc_x is an even number
                    temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                    temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.int8)
                    temp_x1.shape = (lz,ly,lx)
                    temp_x2.shape = (lz,ly,lx)
                    temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]# for ID.* and SignDist.* which have one extra buffer layer around the data
                    temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                    temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                else: # nproc_x is an odd number
                    if nproc_x==1:
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                        temp_x1.shape = (lz,ly,lx)
                        temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                    elif x_axis < nproc_x-2: # i.e. nproc_x = 3 or 5 or 7 ...
                        temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis],dtype=np.int8)
                        temp_x2 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+x_axis+1],dtype=np.int8)
                        temp_x1.shape = (lz,ly,lx)
                        temp_x2.shape = (lz,ly,lx)
                        temp_x1 = temp_x1[1:lz-1,1:ly-1,1:lx-1]
                        temp_x2 = temp_x2[1:lz-1,1:ly-1,1:lx-1]
                        temp_x1 = np.concatenate((temp_x1,temp_x2),axis=2)
                    #end if
                #end if
                if x_axis == 0:
                    ID_x_temp = temp_x1.copy()
                else:
                    ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
                #end if
            #end for
            if nproc_x !=1 and nproc_x%2 != 0:
                temp_x1 = np.fromfile(data_group[z_axis*nproc_y*nproc_x+y_axis*nproc_x+nproc_x-1],dtype=np.int8)
                ID_x_temp = np.concatenate((ID_x_temp,temp_x1),axis=2)
            #end if

            if y_axis == 0:
                ID_y_temp = ID_x_temp.copy()
            else:
                ID_y_temp = np.concatenate((ID_y_temp,ID_x_temp),axis=1)
            #end if
        #end for

        if z_axis == 0:
            ID_z_temp = ID_y_temp.copy()
        else:
            ID_z_temp = np.concatenate((ID_z_temp,ID_y_temp),axis=0)
        #end if
    #end for
    print '**Info: Phase ID data reconstruction completed for '+data_group[0]
    return ID_z_temp.flatten()
#end def


