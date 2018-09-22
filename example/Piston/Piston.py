import numpy

nx=96
ny=24
nz=24
N=nx*ny*nz

mesh=(nx,ny,nz)
data=numpy.ones(mesh,dtype=numpy.int8)

#print(data)
print("Writing piston")
print("Mesh size: "+repr(mesh))

radius = 8
# assign a bubble in the middle
for x in range(0,nx):
  for y in range(0,ny):
    for z in range(0,nz):
       Y = y - ny/2
       Z = z - nz/2
       if Y*Y+Z*Z > radius*radius:
          data[x,y,z]=0
       elif x < 12:
          data[x,y,z]=1
       else:
          data[x,y,z]=2

data.tofile("Piston.raw")
