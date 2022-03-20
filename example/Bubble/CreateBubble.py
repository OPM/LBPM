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
            if (dist < 12.5 ) :
                D[i,j,k] = 2
                
D.tofile("bubble_40x40x40.raw")

