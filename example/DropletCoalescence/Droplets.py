import numpy
import math

nx=400
ny=200
nz=200
N=nx*ny*nz
Radius=64

mesh=(nx,ny,nz)
data=numpy.ones(mesh,dtype=numpy.int8)

#print(data)
print("Create two droplets")
print("Mesh size: "+repr(mesh))
print("Droplet radius: "+repr(Radius))

gap = 6
c1x = nx/2 - gap/2 - Radius
c2x = nx/2 + gap/2 + Radius

# assign a bubble on each side
for x in range(0,200):
  for y in range(0,ny):
    for z in range(0,nz):
       if math.sqrt((x-c1x)*(x-c1x)+(y-ny/2)*(y-ny/2)+(z-nz/2)*(z-nz/2) ) < Radius:
           data[x,y,z]=2

for x in range(200,nx):
  for y in range(0,ny):
    for z in range(0,nz):
       if math.sqrt((x-c2x)*(x-c2x)+(y-ny/2)*(y-ny/2)+(z-nz/2)*(z-nz/2) ) < Radius:
           data[x,y,z]=2

data.tofile("Droplets.raw")
