#!/usr/bin/env python

import os
import numpy as np
import csv

def get_LBM_parameters(csv_file_name,base_name_suffix,image_format,ift=24,Depth=8.8,**kwargs):
    # 'ift': dyne/cm
    # 'Depth': micron
    # Users need to provide the following information
    # 1. surface tension in physical unit
    # 2. physical depth  in physical unit
    # 3. experimental file e.g. *.csv
    # 4. Other optional information, which is summarised in a dictionary
    #    'Para', including: 
    #    Para['D']: characteristic length
    #    Para['alpha']: LBM parameters, controlling the surface tension
    #    Para['fitting']: LBM parameters, extracted from the bubble test


    # Experimental definitions - units are converted in terms of N, m
    micron=1e-6
    cm=1e-2
    dyne=1e-5
    kPa=1e3
    Para={'D':30.0} # characteristic length
    #length=500*micron
    
    # LBM related parameters
    # TODO: for parameters like 'alpha', you should find a way to 
    #       communicate with other functions
    #       It cannot be hard-coded here
    Para['alpha'] = 0.01
    Para['fitting'] = 5.796
    # Check users' input arguments
    for key in kwargs:
        if key in Para.keys():
            Para[key]=kwargs[key]
        else:
            print "Error: "+key+" is not a valid input !\n"
            print "Error: LBM pressure boundaries are not set successfully !"
            return
        #end if
    #end for
    IFT = Para['alpha'] * Para['fitting']

    # 'ift' is the surface tension from the experiment in physical unit
    ift=ift*dyne/cm # This converts the unit of 'ift' to SI unit 
    
    # Process the experimental data
    # 'EXP_data'  : original experimental data
    #               It is a recorded numpy array, which means that its component
    #               can be accessed by 'EXP_data.sw' or 'EXP_data['sw']'
    # About 'EXP_data' see more information from the function 'read_csv'
    EXP_data = read_csv(csv_file_name)
    # A few notes for the experimental data
    # 'Jwn': mean curvature in physical unit (e.g. 1/[micron])
    # 'pwn': pressure difference in physical (e.g. kPa)

    # Overall we need to map the measured physical pressures to get the principle radius of curvature R1 
    # and associated curvature J1, and similarly R2 and J2 in the model's depth direction
    # J1: principal curvature 1 
    # J2: principal curvature 2 along the model's depth direction

    # 'pc' is the dimensionless mean curvature scaled by the length scale 'D32'
    # 'pc' is extracted from experimentally measured mean curvature 'Jwn'
    pc = Para['D']*micron*EXP_data.Jwn*(1.0/micron)

    # Alternatively, the dimensionless mean curvature can also be extracted from the 
    # experimentally measured pressure difference (i.e. capillary pressure)
    # 'pwn' is the dimensionless mean curvature scaled by the length scale 'D32' 
    # pwn = Para['D']*micron*EXP_data.pwn*kPa/ift 

    # Curvature is fixed in the micromodel "depth" direction
    J2 = Para['D']*micron/(Depth*micron/2.0) 

    # infer the curvature in the other direction
    J1 = pc-J2

    # determine the LBM pressure difference
    dp=(J1+J2)*IFT/Para['D']
    # determine the boundary pressure values for Color.in
    din  = 1.0+0.5*dp 
    dout = 1.0-0.5*dp
    
    # Generate names of segmented image data for Segmented.in file
    # Need the input 'base_name_suffix' (e.g. '_segmented.raw')
    data_name_for_Segmented_in = EXP_data.image_name.copy()
    data_name_for_Segmented_in = data_name_for_Segmented_in.replace(image_format,base_name_suffix)

    return (data_name_for_Segmented_in,din,dout)

#end def

def write_Domain_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'nprocx':1}
    Para['nprocy']=2
    Para['nprocz']=2
    Para['nx']=1
    Para['ny']=2
    Para['nz']=2
    Para['nspheres']=0 # deprecated
    Para['Lx']=10
    Para['Ly']=500
    Para['Lz']=500

    for key in kwargs:
        if key in Para.keys():
            Para[key]=kwargs[key]
        else:
            print "Error: "+key+" is not a valid input !\n"
            print "Error: Domain.in file is not writen."
            return
        #end if
    #end for
    
    # Check if the dompositoin requirement is satisfied
    if (Para['nx']*Para['nprocx']>Para['Lx']) or \
       (Para['ny']*Para['nprocy']>Para['Ly']) or \
       (Para['nz']*Para['nprocz']>Para['Lz']):
            print "Error: The decomposed size does not match the total domain size!"
            return
    #end if

    # write the Domain.in file
    ParamFile = open("Domain.in","w")
    ParamFile.write("%i "   % Para['nprocx'])
    ParamFile.write("%i "   % Para['nprocy'])
    ParamFile.write("%i\n" % Para['nprocz'])
    ParamFile.write("%i "   % Para['nx'])
    ParamFile.write("%i "   % Para['ny'])
    ParamFile.write("%i\n" % Para['nz'])
    ParamFile.write("%i\n" % Para['nspheres'])
    ParamFile.write("%.1f "   % Para['Lx'])
    ParamFile.write("%.1f "   % Para['Ly'])
    ParamFile.write("%.1f\n" % Para['Lz'])
    ParamFile.close()
    print "A Domain.in file is created."
#end def


def write_Segment_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'Seg_data_name':'Micromodel_1_segmented.raw'}
    Para['Nx']=10
    Para['Ny']=500
    Para['Nz']=500
    Para['xStart']=0
    Para['yStart']=0
    Para['zStart']=0

    for key in kwargs:
        if key in Para.keys():
            Para[key]=kwargs[key]
        else:
            print "Error: "+key+" is not a valid input !\n"
            print "Error: Segmented.in file is not writen."
            return
        #end if
    #end for

    ParamFile = open("Segmented.in","w")
    ParamFile.write("%s\n" % Para['Seg_data_name'])
    ParamFile.write("%i "  % Para['Nx'])
    ParamFile.write("%i "  % Para['Ny'])
    ParamFile.write("%i\n" % Para['Nz'])
    ParamFile.write("%i "  % Para['xStart'])
    ParamFile.write("%i "  % Para['yStart'])
    ParamFile.write("%i\n" % Para['zStart'])
    ParamFile.close()
    print "A Segmented.in file is created."
#end def

def write_Color_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'tau':0.7}
    Para['alpha']=0.01
    Para['beta']=0.95
    Para['phi_solid']=-1.0
    Para['saturation']=0.0
    Para['Fx']=0.0
    Para['Fy']=0.0
    Para['Fz']=0.0
    Para['Restart']=0
    Para['pBC']=1
    Para['din']=1.001
    Para['dout']=0.999
    Para['maxtime']=100005
    Para['interval']=2000
    Para['tolerance']=1e-5

    for key in kwargs:
        if key in Para.keys():
            Para[key]=kwargs[key]
        else:
            print "Error: "+key+" is not a valid input !\n"
            print "Error: Color.in file is not writen."
            return
        #end if
    #end for
    
    # write the color.in file
    ParamFile = open("Color.in","w")
    ParamFile.write("%f\n" % Para['tau'])
    ParamFile.write("%f " % Para['alpha'])
    ParamFile.write("%f " % Para['beta'])
    ParamFile.write("%f\n" % Para['phi_solid'])
    ParamFile.write("%f\n" % Para['saturation'])
    ParamFile.write("%f " % Para['Fx'])
    ParamFile.write("%f " % Para['Fy'])
    ParamFile.write("%f\n" % Para['Fz'])
    ParamFile.write("%i " % Para['Restart'])
    ParamFile.write("%i " % Para['pBC'])
    ParamFile.write("%f " % Para['din'])
    ParamFile.write("%f\n" % Para['dout'])
    ParamFile.write("%i " % Para['maxtime'])
    ParamFile.write("%i " % Para['interval'])
    ParamFile.write("%e\n" % Para['tolerance'])
    ParamFile.close()
    print "A Color.in file is created."
#end def

def read_csv(csv_file_name):
    #TODO Haven't thought about a better way of doing this
    # Right now I just hard-code the possible variables from the experiment
    # which means users are forced to prepare their *.csv file as listed here
    image_name = []
    sw = []
    Jwn = []
    pwn = []
    with open(csv_file_name,"r") as f:
        for line in f:
            reader=csv.reader(f,delimiter=' ') #Note the header is skipped by default 
            for row in reader:
                image_name.append(row[0])
                sw.append(float(row[1]))
                Jwn.append(float(row[2]))
                pwn.append(float(row[3]))
            #end for
        #end for
    #end with
    
    # Convter the list to numpy array
    image_name_array = np.asarray(image_name,dtype=str)
    sw_array  = np.asarray(sw,dtype=np.float32)
    Jwn_array = np.asarray(Jwn,dtype=np.float32)
    pwn_array = np.asarray(pwn,dtype=np.float32)

    # Restore the original shape of the experimental data
    experiment_data = np.column_stack((image_name_array,sw_array,Jwn_array,pwn_array))
    # Unfortunately the column stack will cast all float data to strings
    experiment_data = np.core.records.fromrecords(experiment_data,names='image_name,sw,Jwn,pwn')
    dt=experiment_data.dtype.descr
    dt[1] = (dt[1][0],'float32')
    dt[2] = (dt[2][0],'float32')
    dt[3] = (dt[3][0],'float32')
    experiment_data = experiment_data.astype(dt)
    return experiment_data
#end def






