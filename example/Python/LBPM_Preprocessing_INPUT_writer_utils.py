#!/usr/bin/env python

import os
import numpy as np

def write_Domain_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'nprocx':1}
    Para['nprocy']=1
    Para['nprocz']=1
    Para['nx']=1
    Para['ny']=1
    Para['nz']=1
    Para['nspheres']=0 # deprecated
    Para['Lx']=100
    Para['Ly']=100
    Para['Lz']=100

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
    print "**Info: A Domain.in file is created."
#end def


def write_Segment_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'Seg_data_name':'Micromodel_1_segmented.raw'}
    Para['Nx']=100
    Para['Ny']=100
    Para['Nz']=100
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
    print "**Info: A Segmented.in file is created."
#end def

def write_Color_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'tau':0.7}
    Para['alpha']=0.001
    Para['beta']=0.95
    Para['phi_solid']=-0.98
    Para['saturation']=0.7 # wetting phase saturation, depricated
    Para['Fx']=0.0
    Para['Fy']=0.0
    Para['Fz']=0.0
    Para['InitialCondition']=0
    Para['BoundaryCondition']=0
    Para['din']=0.0
    Para['dout']=0.0
    Para['TimeStepMax']=100005
    Para['Restart_interval']=2000
    Para['tolerance']=1e-5
    Para['Blob_analysis_interval'] = 1000 #NOTE: not used

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
    ParamFile.write("%.3f\n" % Para['tau'])
    ParamFile.write("%.3e " % Para['alpha'])
    ParamFile.write("%.3f " % Para['beta'])
    ParamFile.write("%.3f\n" % Para['phi_solid'])
    ParamFile.write("%.2f\n" % Para['saturation'])
    ParamFile.write("%.3e " % Para['Fx'])
    ParamFile.write("%.3e " % Para['Fy'])
    ParamFile.write("%.3e\n" % Para['Fz'])
    ParamFile.write("%i " % Para['Restart'])
    ParamFile.write("%i " % Para['pBC'])
    ParamFile.write("%.3f " % Para['din'])
    ParamFile.write("%.3f\n" % Para['dout'])
    ParamFile.write("%i " % Para['maxtime'])
    ParamFile.write("%i " % Para['interval'])
    ParamFile.write("%.2e\n" % Para['tolerance'])
    ParamFile.close()
    print "**Info: A Color.in file is created."
#end def

def write_Color_Macro_in_file(**kwargs):
    # For most of the parameters in the Color.in file
    # you wouldn't need to worry about them as by default 
    # they are all hard-coded here 
    Para={'tau1':1.0}
    Para['tau2'] = 1.0
    Para['alpha']=0.001
    Para['beta']=0.95
    Para['phi_solid']=-0.98
    Para['saturation']=0.7 # wetting phase saturation, depricated
    Para['Fx']=0.0
    Para['Fy']=0.0
    Para['Fz']=0.0
    Para['InitialCondition']=0
    Para['BoundaryCondition']=0
    Para['din']=0.0
    Para['dout']=0.0
    Para['TimeStepMax']=100005
    Para['Restart_interval']=2000
    Para['tolerance']=1e-5
    Para['Blob_analysis_interval'] = 1000 #NOTE: not used

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
    ParamFile.write("%.3f " % Para['tau1'])
    ParamFile.write("%.3f " % Para['tau2'])
    ParamFile.write("%.3e " % Para['alpha'])
    ParamFile.write("%.3f " % Para['beta'])
    ParamFile.write("%.3f\n" % Para['phi_solid'])
    ParamFile.write("%.2f\n" % Para['saturation'])
    ParamFile.write("%.3e " % Para['Fx'])
    ParamFile.write("%.3e " % Para['Fy'])
    ParamFile.write("%.3e\n" % Para['Fz'])
    ParamFile.write("%i " % Para['InitialCondition'])
    ParamFile.write("%i " % Para['BoundaryCondition'])
    ParamFile.write("%.3f " % Para['din'])
    ParamFile.write("%.3f\n" % Para['dout'])
    ParamFile.write("%i " % Para['TimeStepMax'])
    ParamFile.write("%i " % Para['Restart_interval'])
    ParamFile.write("%.2e\n" % Para['tolerance'])
    #ParamFile.write("%i\n" % Para['Blob_analysis_interval'])
    ParamFile.close()
    print "**Info: A Color.in file is created."
#end def

