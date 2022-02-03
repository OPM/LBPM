import numpy as np
import matplotlib.pylab as plt

D=np.ones((80,80,40),dtype="uint8")

cx = 40
cy = 40
cz = 0

for i in range(0,80):
    for j in range (0, 80):
        for k in range (0,40):
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < 25) :
                D[i,j,k] = 2
            if k < 4 : 
                D[i,j,k] = 0
                
D.tofile("droplet_40x80x80.raw")

