import numpy as np

N = 1024

data = np.random.randint(low=1,high=3,size=(N,N,N),dtype=np.uint8)

data.tofile("dense_1024x1024x1024.raw")


