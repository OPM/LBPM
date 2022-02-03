*************************
Measuring Contact Angles
*************************

LBPM includes specialized data analysis capabilities for two-fluid systems. While these
components are generally designed for in situ analysis of simulation data, they can also
be applied independently to analyze 3D image data. In this example we consider applying
the analysis tools implemented in ``lbpm_TwoPhase_analysis``, which are designed to
analyze two-fluid configurations in porous media. The numerical implementation used to
construct the common line are described in ( https://doi.org/10.1016/j.advwatres.2006.06.010 ).
Methods used to measure the contact angle are described in ( https://doi.org/10.1017/jfm.2016.212 ).

Source files for the example are included in the LBPM repository
in the directory ``examples/Droplet``. A simple python code is included
to set up a fluid droplet on a flat surface

.. code:: python
	  
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


The input file provided below will specify how the analysis should be performed. The name and the dimensions of the input
file are provided in the ``Domain`` section, as with other LBPM simulators. For large images, additional processors can be
used to speed up the analysis or take advantage of distributed memory. The ``ReadValues`` list should be used to specify the
labels to use for analysis. The first label will be taken to be the solid.  The second label will be taken to be the fluid
to analyze, which in this case will be the droplet labeled with ``2`` above. 
	  
.. code:: bash
	  
       Domain {
	  Filename = "droplet_40x80x80.raw"
	  nproc = 1, 1, 1     // Number of processors (Npx,Npy,Npz)
	  n = 40, 80, 80      // Size of local domain (Nx,Ny,Nz)
	  N = 40, 80, 80      // size of the input image
	  voxel_length = 1.0 
	  BC = 0              // Boundary condition type
	  ReadType = "8bit"
	  ReadValues = 0, 2, 1
	  WriteValues = 0, 2, 1
      }
      Visualization {
      }


The analysis can be launched as ``mpirun -np 1 $LBPM_DIR/lbpm_TwoPhase_analysis input.db``. Output should appear as
follows:

.. code:: bash
 
	   Input data file: input.db
	   voxel length = 1.000000 micron 
	   Input media: droplet_40x80x80.raw
	   Relabeling 3 values
	   oldvalue=0, newvalue =0 
	   oldvalue=2, newvalue =2 
	   oldvalue=1, newvalue =1 
	   Dimensions of segmented image: 40 x 80 x 80 
	   Reading 8-bit input data 
	   Read segmented data from droplet_40x80x80.raw 
	   Label=0, Count=25600 
	   Label=2, Count=25773 
	   Label=1, Count=204627 
	   Distributing subdomains across 1 processors 
	   Process grid: 1 x 1 x 1 
	   Subdomain size: 40 x 80 x 80 
	   Size of transition region: 0 
	   Media porosity = 0.900000 
	   Initialized solid phase -- Converting to Signed Distance function 
	   Initialized fluid phase -- Converting to Signed Distance function 
	   Computing Minkowski functionals 

The ``TwoPhase`` analysis class will generate signed distance functions for the solid and fluid surfaces.
Using the distance functions, the interfaces and common line will be constructed. Contact angles are logged
for each processor, e.g. ``ContactAngle.00000.csv``, which specifies the x, y, z coordinates for each measurement
along with the cosine of the contact angle. Averaged measures (determined over the entire input image)
are logged to ``geometry.csv``

* ``sw`` -- water saturation
* ``awn`` -- surface area of meniscus between wn fluids
* ``ans`` -- surface area between fluid n and solid
* ``aws`` -- surface area between fluid w and solid
* ``Jwn`` -- integral of mean curvature of meniscus
* ``Kwn`` -- integral of Gaussian curvature of meniscus
* ``lwns`` -- length of common line
* ``cwns`` -- average contact angle
* ``KGws`` -- geodesic curvature of common line relative to ws surface
* ``KGwn`` -- geodesic curvature of common line relative to wn surface
* ``Gwnxx`` -- orientation tensor component for wn surface
* ``Gwnyy`` -- orientation tensor component for wn surface
* ``Gwnzz`` -- orientation tensor component for wn surface
* ``Gwnxy`` -- orientation tensor component for wn surface
* ``Gwnxz`` -- orientation tensor component for wn surface
* ``Gwnyz`` -- orientation tensor component for wn surface
* ``Gwsxx`` -- orientation tensor component for ws surface
* ``Gwsyy`` -- orientation tensor component for ws surface
* ``Gwszz`` -- orientation tensor component for ws surface
* ``Gwsxy`` -- orientation tensor component for ws surface
* ``Gwsxz`` -- orientation tensor component for ws surface
* ``Gwsyz`` -- orientation tensor component for ws surface
* ``Gnsxx`` -- orientation tensor component for ns surface
* ``Gnsyy`` -- orientation tensor component for ns surface
* ``Gnszz`` -- orientation tensor component for ns surface
* ``Gnsxy`` -- orientation tensor component for ns surface
* ``Gnsxz`` -- orientation tensor component for ns surface
* ``Gnsyz`` -- orientation tensor component for ns surface
* ``trawn`` -- trimmed surface area for meniscus (one voxel from solid)
* ``trJwn`` -- mean curvature for trimmed meniscus
* ``trRwn`` -- radius of curvature for trimmed meniscus
* ``Vw`` -- volume of fluid w
* ``Aw`` -- boundary surface area for fluid w
* ``Jw`` -- integral of mean curvature for fluid w
* ``Xw`` -- Euler characteristic for fluid w
* ``Vn`` -- volume of fluid n
* ``An``  -- boundary surface area for fluid n
* ``Jn`` -- integral of mean curvature for fluid n
* ``Xn`` -- Euler characteristic for fluid n


  
