###############################################################################
MRT model
###############################################################################

The LBPM single fluid model is implemented by combining a multi-relaxation time (MRT) D3Q19
lattice Boltzmann equation (LBE) to solve for the momentum transport, recovering the Navier-Stokes
equations to second order based on the Chapman-Enskog expansion. The MRT model is used to assess the
permeability of digital rock images in either the Darcy or non-Darcy flow regimes. 

A typical command to launch the LBPM color simulator is as follows

```
mpirun -np $NUMPROCS lbpm_permeability_simulator input.db
```

where ``$NUMPROCS`` is the number of MPI processors to be used and ``input.db`` is
the name of the input database that provides the simulation parameters.
Note that the specific syntax to launch MPI tasks may vary depending on your system.
For additional details please refer to your local system documentation.

***************************
Model parameters
***************************

The essential model parameters for the single-phase MRT model are

- ``tau`` -- control the fluid viscosity -- :math:`0.7 < \tau < 1.5`

The kinematic viscosity is given by

.. math::
   :nowrap:

     $$
       \nu = \frac{1}{3} \Big( \tau - \frac 12 \Big)
     $$

****************************
Model Formulation
****************************

The LBE governing momentum transport is defined based on a MRT relaxation based on the D3Q19 discrete
velocity set, which determines the values :math:`\bm{\xi}_q`

.. math::
   :nowrap:

   $$
      f_q(\bm{x}_i + \bm{\xi}_q \delta t,t + \delta t) - f_q(\bm{x}_i,t) = \sum^{Q-1}_{k=0} M^{-1}_{qk} \lambda_{k} (m_k^{eq}-m_k) + w_q \bm{\xi}_q \cdot \frac{\bm{F}}{c_s^2} \;,
   $$

Where :math:`\bm{F}` is an external body force and :math:`c_s^2 = 1/3` is the speed of sound for the LB model.
The moments are linearly indepdendent functions of the distributions:

.. math::
   :nowrap:

   $$
      m_k = \sum_{q=0}^{18} M_{qk} f_q\;.
   $$


The non-zero equilibrium moments are

.. math::
   :nowrap:

   $$
     m_1^{eq} = 19\frac{j_x^2+j_y^2+j_z^2}{\rho} - 11\rho \;,
   $$     

.. math::
   :nowrap:

   $$
     m_2^{eq} = 3\rho - \frac{11}{2} \frac{j_x^2+j_y^2+j_z^2}{\rho} \;,
   $$     

.. math::
   :nowrap:


   $$
     m_4^{eq} = -\frac 2 3 j_x \;,
   $$
   
.. math::
   :nowrap:

   $$
     m_6^{eq} = -\frac 2 3 j_y \;,
   $$

.. math::
   :nowrap:

   $$
     m_8^{eq} = -\frac 2 3 j_z \;,
   $$

.. math::
   :nowrap:

   $$     
     m_9^{eq} = \frac{2j_x^2-j_y^2-j_z^2}{\rho}\;,
   $$     

.. math::
   :nowrap:

   $$
     m_{10}^{eq} = -\frac{2j_x^2-j_y^2-j_z^2)}{2\rho} \;,
   $$     
   
.. math::
   :nowrap:

   $$     
     m_{11}^{eq} = \frac{j_y^2-j_z^2}{\rho} \;, 
   $$     

.. math::
   :nowrap:

   $$
     m_{12}^{eq} = -\frac{j_y^2-j_z^2}{2\rho} \;,
   $$     

   
.. math::
   :nowrap:

   $$     
     m_{13}^{eq} = \frac{j_x j_y}{\rho} \;, 
   $$     

.. math::
   :nowrap:

   $$     
     m_{14}^{eq} = \frac{j_y j_z}{\rho} \;, 
   $$     

.. math::
   :nowrap:

   $$     
     m_{15}^{eq} = \frac{j_x j_z}{\rho} \;, 
   $$

The relaxation parameters are determined based on the relaxation time :math:`\tau`

.. math::
   :nowrap:

   $$
     \lambda_1 =  \lambda_2=  \lambda_9 = \lambda_{10}= \lambda_{11}= \lambda_{12}= \lambda_{13}= \lambda_{14}= \lambda_{15} = s_\nu = \frac{1}{\tau} \;,
   $$
   
.. math::
   :nowrap:
      
    $$
     \lambda_{4}= \lambda_{6}= \lambda_{8} = \lambda_{16} = \lambda_{17} = \lambda_{18}= \frac{8(2-s_\nu)}{8-s_\nu} \;,
   $$


  
****************************
Boundary Conditions
****************************

The following external boundary conditions are supported by ``lbpm_permeability_simulator``
and can be set by setting the ``BC`` key values in the ``Domain`` section of the
input file database

- ``BC = 0`` -- fully periodic boundary conditions
- ``BC = 3`` -- constant pressure boundary condition
- ``BC = 4`` -- constant volumetric flux boundary condition

For ``BC = 0`` any mass that exits on one side of the domain will re-enter at the other
side. If the pore-structure for the image is tight, the mismatch between the inlet and
outlet can artificially reduce the permeability of the sample due to the blockage of
flow pathways at the boundary. LBPM includes an internal utility that will reduce the impact
of the boundary mismatch by eroding the solid labels within the inlet and outlet layers
(https://doi.org/10.1007/s10596-020-10028-9) to create a mixing layer.
The number mixing layers to use can be set using the key values in the ``Domain`` section
of the input database

- ``InletLayers  = 5`` -- set the number of mixing layers to ``5`` voxels at the inlet
- ``OUtletLayers  = 5`` -- set the number of mixing layers to ``5`` voxels at the outlet

For the other boundary conditions a thin reservoir of fluid  (default ``3`` voxels)
is established at either side of the domain. The inlet is defined as the boundary face
where ``z = 0`` and the outlet is the boundary face where ``z = nprocz*nz``. By default a
reservoir of fluid A is established at the inlet and a reservoir of fluid B is established at
the outlet, each with a default thickness of three voxels. To over-ride the default label at
the inlet or outlet, the ``Domain`` section of the database may specify the following key values

- ``InletLayerPhase = 2`` -- establish a reservoir of component B at the inlet
- ``OutletLayerPhase = 1`` -- establish a reservoir of component A at the outlet

  ****************************
Example Input File
****************************

.. code-block:: c

   MRT {
      tau = 1.0
      F = 0.0, 0.0, 1.0e-5
      timestepMax = 2000
      tolerance = 0.01
   }
   Domain {
      Filename = "Bentheimer_LB_sim_intermediate_oil_wet_Sw_0p37.raw"  
      ReadType = "8bit"      // data type
      N = 900, 900, 1600     // size of original image
      nproc = 2, 2, 2        // process grid
      n = 200, 200, 200      // sub-domain size
      offset = 300, 300, 300 // offset to read sub-domain
      voxel_length = 1.66    // voxel length (in microns)
      ReadValues = 0, 1, 2   // labels within the original image
      WriteValues = 0, 1, 2  // associated labels to be used by LBPM
      InletLayers = 0, 0, 10 // specify 10 layers along the z-inlet
      BC = 0                 // boundary condition type (0 for periodic)
   }
   Visualization {
   }
