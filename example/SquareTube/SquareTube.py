import numpy

nx=96
ny=24
nz=24
N=nx*ny*nz

mesh=(nx,ny,nz)
data=numpy.ones(mesh,dtype=numpy.int8)

LabelTop=-1
LabelBottom=-2
LabelLeft=-3
LabelRight=-4

#print(data)
print("Writing parllel plates")
print("Mesh size: "+repr(mesh))
print("Top label: "+repr(LabelTop))
print("Bottom label: "+repr(LabelBottom))
print("Left label: "+repr(LabelLeft))
print("Right label: "+repr(LabelRight))

# assign a bubble in the middle
for x in range(12,72):
  for y in range(0,ny):
    for z in range(0,nz):
       data[x,y,z]=2

# solid walls with different component labels
for x in range(0,nx):
  for y in range(0,ny):
    data[x,y,0]=LabelBottom
    data[x,y,nz-1]=LabelTop


for x in range(0,nx):
  for z in range(0,nz):
    data[x,0,z]=LabelLeft
    data[x,ny-1,z]=LabelRight


data.tofile("SquareTube.raw")
