import numpy as np
import matplotlib.pylab as plt

D=np.ones((40,40,40),dtype="uint8")

cx = 20
cy = 20
cz = 20
radius = 8

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < radius ) :
                D[i,j,k] = 2

D.tofile("cell_40x40x40.raw")

                
C1=np.zeros((40,40,40),dtype="double")
C2=np.zeros((40,40,40),dtype="double")

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            #outside the cell
            C1[i,j,k] = 125.0e-6  # Na
            C2[i,j,k] = 125.0e-6  # Cl
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            # inside the cell
            if (dist < radius ) :
                C1[i,j,k] = 5.0e-6
                C2[i,j,k] = 5.0e-6

C1.tofile("cell_concentration_Na_40x40x40.raw")
C2.tofile("cell_concentration_Cl_40x40x40.raw")

