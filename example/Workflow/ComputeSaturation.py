import sys
import numpy as np
import matplotlib.pylab as plt

FILENAME=sys.argv[1]
Nx=int(sys.argv[2])
Ny=int(sys.argv[3])
Nz=int(sys.argv[4])

# read the input image
Output = np.fromfile(FILENAME,dtype = np.uint8)
Output.shape = (Nz,Ny,Nx)

Oil=np.count_nonzero(Output==1)
Water=np.count_nonzero(Output==2)
Sw=Water/(Oil+Water)

Porosity=1.0-(Oil+Water)/(Nx*Ny*Nz)

print(FILENAME,"Porosity=", Porosity)

SaturationProfile=np.zeros(Nz)
PorosityProfile=np.zeros(Nz)
# Compute saturation slice by slice 
for idx in range(0, Nz):
   Slice = Output[idx,:,:]
   Oil=np.count_nonzero(Slice==1)
   Water=np.count_nonzero(Slice==2)
   SaturationProfile[idx]=Water/(Oil+Water)
   PorosityProfile[idx]=(Oil+Water)/(Nx*Ny)
   

plt.figure()
plt.plot(SaturationProfile)
plt.xlabel('Position (z)')
plt.ylabel('Water Saturation')
plt.show()
