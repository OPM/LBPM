import numpy as np
import matplotlib.pylab as plt

D=np.ones((40,40,40),dtype="uint8")

cx = 20
cy = 20
cz = 20

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < 15.5 ) :
                D[i,j,k] = 2

D.tofile("cell_40x40x40.raw")

                
C1=np.zeros((40,40,40),dtype="double")
C2=np.zeros((40,40,40),dtype="double")

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            C1[i,j,k] = 4.0e-6
            C2[i,j,k] = 150.0e-6
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < 15.5 ) :
                C1[i,j,k] = 140.0e-6
                C2[i,j,k] = 10.0e-6

                
C1.tofile("cell_concentration_K_40x40x40.raw")
C2.tofile("cell_concentration_Na_40x40x40.raw")

