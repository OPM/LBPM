*****************
Input Domain
*****************

LBPM provides a flexible framework to ingest 3D image data.
To illustrate the basic capabilities, this tutorial considers a quasi-2D
flow cell. Source files for the example are included in the LBPM repository
in the directory ``examples/DiscPack``. A simple python code is included
to set up the flow domain.

Based on LBPM convention, external boundary conditions are applied in the
z-direction. This means that the domain should be set up so that the direction
you want to set boundary conditions is aligned with the z-axis. For the quasi-2D
example, a depth of ``3`` voxels is used for the x-direction. *Based on LBPM
internal data structures at least three voxels must be provided in each direction*
The specified domain decomposition must also observe this rule.

Image data is stored internally within LBPM as signed 8-bit binary data. This means that
up to 256 labels can be provided in the input image. LBPM convention takes all
non-positive labels to be immobile (treated as solid). In this example, the solid regions
are assigned a value of ``0``.  It is possible to provide up to ``128`` different labels
for the solid. Also, note that python supports only the unsigned 8-bit datatype. For the unsigned data
type, labels assigned values ``128,...255`` in python will correspond to labels
``-127,...-1`` when read in as type ``signed char`` within LBPM. 

.. code:: python

	  import numpy as np
	  import matplotlib.pylab as plt
	  import pandas as pd
	  # Set the size of the domain
	  Nx=3
	  Ny=128
	  Nz=128
	  D=pd.read_csv("discs.csv",sep=" ")
	  ID = np.ones(Nx*Ny*Nz,dtype='uint8')
	  ID.shape = (Nz,Ny,Nx)
	  # Set the solid labels
	  for idx in range(len(D)):
             cx=D['cx'][idx] / dx
             cy=D['cy'][idx] /dx
             r=D['r'][idx] /dx
             for i in range(0,Nz):
                for j in range(0,Ny):
                   if ( (cx-i)*(cx-i) + (cy-j)*(cy-j) < r*r ):
                      ID[i,j,0] = 0
                      ID[i,j,1] = 0
                      ID[i,j,2] = 0
        # write input file to disc
	ID.tofile("discs_3x128x128.raw")


