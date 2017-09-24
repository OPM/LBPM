#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt


# Check if there is a proper command line argument
if len(sys.argv) !=2:
    sys.stderr.write('**Error: Usage: ' + sys.argv[0] + ' <Color.in>\n')
    sys.exit()
# end if

## *********************** Read 'Color.in' ***************************** ##
print "**Info: Reading Color.in file."
f = open(sys.argv[1],'r')
lines = f.readlines()
tau=float(lines[0])
alpha, beta, phi_s = np.fromstring(lines[1].splitlines()[0],dtype=np.float64, sep=' ')
# Read the 3-D body forces
Fx, Fy, Fz = np.fromstring(lines[3].splitlines()[0],dtype=np.float64, sep=' ')
# There can be more to be read from Color.in file
f.close()
## ******************** END: Read 'Color.in' *************************** ##

# Load the data file
print "**Info: Reading timelog.tcat file."
data_all = np.genfromtxt('timelog.tcat',names=True)
# Load individual data sets
pw_all=data_all['pw']
#pw=pw_all[np.isfinite(pw_all)]
sw=data_all['sw']
sw=sw[np.isfinite(pw_all)]
vanz=data_all['vanz']
vanz=vanz[np.isfinite(pw_all)]
time_step=data_all['time']
time_step=time_step[np.isfinite(pw_all)]
Jwn = data_all['Jwn']
Jwn=Jwn[np.isfinite(pw_all)]
## ********************************** LBPM-WIA related parameters *********************************** ##
# (lattice) interfacial tension
#TODO: need to load Color.in file
IFT_conv = 5.796
IFT = alpha*IFT_conv
viscosity = (tau-0.5)/3.0 # lattice dynamic viscosity
cos_theta=phi_s #NOTE: this is only an approximation

# Plot Capillary number
#Ca=vanz*viscosity/IFT/cos_theta
Ca=vanz*viscosity/IFT
plt.figure(1)
plt.subplot(1,2,1)
plt.plot(time_step,Ca,'ks--',label='Total Ca')
plt.plot(np.array([0,time_step.max()]),np.array([0,0]),'r-',linewidth=2)
#plt.plot(time_step,vanz,'ro--')
plt.grid(True)
plt.legend(loc='best')
plt.xlabel('Time step (LBM)')
plt.ylabel('Capillary number (w.r.t non-wetting phase)')

plt.subplot(1,2,2)
plt.semilogy(time_step[Ca>0],Ca[Ca>0],'b>--',label='Positive Ca')
plt.semilogy(time_step[Ca<0],np.abs(Ca[Ca<0]),'ro--',label='Negative Ca')
#plt.plot(time_step,vanz,'ro--')
plt.grid(True)
plt.legend(loc='best')
plt.xlabel('Time step (LBM)')
plt.ylabel('Capillary number (w.r.t non-wetting phase)')

plt.figure(2)
plt.plot(time_step,sw,'ro--')
plt.grid(True)
plt.xlabel('Time step (LBM)')
plt.ylabel('Wetting phase saturation Sw')


plt.figure(3)
plt.plot(time_step,-1*Jwn,'b>--')
plt.grid(True)
plt.xlabel('Time step (LBM)')
plt.ylabel('Interfacial mean curvature (-1.0*Jwn)')


plt.show()











