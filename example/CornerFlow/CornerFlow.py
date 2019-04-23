import numpy
import math

nx=128
ny=64
nz=64
N=nx*ny*nz

mesh=(nx,ny,nz)
data=numpy.ones(mesh,dtype=numpy.int8)

#ylinder centers
c1x=-1.0
c1y=0.0
c2x=1.0
c2y=0.0
c3x=0.0
c3y=math.sqrt(3.0)

# domain size
x0=-0.5
y0=0.0
L=1.0
#2.0*(1.0-1.0/math.sqrt(2.0))
dx=L/nz
radius=1.0

print("Mesh size: "+repr(mesh))
print("Length: "+repr(L))
print("dx: "+repr(dx))
print("box corner: "+repr(x0)+","+repr(y0))
print("cylinder 1: "+repr(c1x)+","+repr(c1y))
print("cylinder 2: "+repr(c2x)+","+repr(c2y))
print("cylinder 3: "+repr(c3x)+","+repr(c3y))

#print(data)
count=0
# set up the tube with corner flow
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
       z = i*dx
       x = j*dx+x0
       y = k*dx+y0
       label=2
       d1 = math.sqrt((x-c1x)*(x-c1x)+(y-c1y)*(y-c1y)) - radius
       d2 = math.sqrt((x-c2x)*(x-c2x)+(y-c2y)*(y-c2y)) - radius
       d3 = math.sqrt((x-c3x)*(x-c3x)+(y-c3y)*(y-c3y)) - radius
       if d1 < 0.0:
         label=0
       if d2 < 0.0:
         label=0
       if d3 < 0.0:
         label=0
       data[i,j,k]=label
       if label==2:
         count+=1

print(count)

data.tofile("CornerFlow.raw")
