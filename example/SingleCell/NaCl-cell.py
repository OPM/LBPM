import numpy as np
import matplotlib.pylab as plt

Nx = 64
Ny = 64
Nz = 64
cx = Nx/2
cy = Ny/2
cz = Nz/2
radius = 12

D=np.ones((Nx,Ny,Nz),dtype="uint8")

for i in range(0, Nx):
    for j in range (0, Ny):
        for k in range (0,Nz):
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < radius ) :
                D[i,j,k] = 2

D.tofile("cell_64x64x64.raw")

                
C1=np.zeros((Nx,Ny,Nz),dtype="double")
C2=np.zeros((Nx,Ny,Nz),dtype="double")

for i in range(0, Nx):
    for j in range (0, Ny):
        for k in range (0,Nz):
            #outside the cell
            C1[i,j,k] = 125.0e-6  # Na
            C2[i,j,k] = 125.0e-6  # Cl
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            # inside the cell
            if (dist < radius ) :
                C1[i,j,k] = 110.0e-6
                C2[i,j,k] = 110.0e-6

C1.tofile("cell_concentration_Na_64x64x64.raw")
C2.tofile("cell_concentration_Cl_64x64x64.raw")

