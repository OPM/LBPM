###############################################################################
Color model
###############################################################################

The LBPM color model is implemented by combining a multi-relaxation time D3Q19
lattice Boltzmann equation (LBE) to solve for the momentum transport with two D3Q7
LBEs for the mass transport. The color model will obey strict mass and momentum
conservation while minimizing diffusive fluxes across the interface between fluids.
The color model is a good choice for modeling dense fluids that are strongly immiscible
(e.g. water-oil systems). Due to the strong anti-diffusion in the interface region,
the color model is not suitable for modeling processes such as Ostwald ripening that
depend on diffusive fluxes between fluid phases.

A typical command to launch the LBPM color simulator is as follows

```
mpirun -np $NUMPROCS lbpm_color_simulator input.db
```

where ``$NUMPROCS`` is the number of MPI processors to be used and ``input.db`` is
the name of the input database that provides the simulation parameters.
Note that the specific syntax to launch MPI tasks may vary depending on your system.
For additional details please refer to your local system documentation.

****************************
Simulation protocols
****************************

Simulation protocols are designed to make it simpler to design and execute common
computational experiments. Protocols will automatically determine boundary conditions
needed to perform a particular simulation. LBPM will internall set default simulation paramaters
that can be over-ridden to develop customized simulations.

.. toctree::
   :glob:
   :maxdepth: 2

   protocols/*

****************************
Analysis capabilities
****************************
   
.. toctree::
   :glob:
   :maxdepth: 2

   analysis/*


***************************
Model parameters
***************************

The essential model parameters for the color model are

- ``alpha`` -- control the interfacial tension between fluids -- :math:`0 < \alpha < 0.01`
- ``beta`` -- control the width of the interface -- :math:`\beta < 1`
- ``tauA`` -- control the viscosity of fluid A -- :math:`0.7 < \tau_A < 1.5`
- ``tauB`` -- control the viscosity of fluid B -- :math:`0.7 < \tau_B < 1.5`
- ``rhoA`` -- control the viscosity of fluid A -- :math:`0.05 < \rho_A < 1.0`
- ``rhoB`` -- control the viscosity of fluid B -- :math:`0.05 < \rho_B < 1.0`

****************************
Model Formulation
****************************


The relaxation parameters are determined from the relaxation time:
.. math::
   :nowrap:


   $$
   \\begin{eqnarray}
     \lambda_1 =  \lambda_2=  \lambda_9 = \lambda_{10}= \lambda_{11}= \lambda_{12}= \lambda_{13}= \lambda_{14}= \lambda_{15} = s_\nu\;, \\
     \lambda_{4}= \lambda_{6}= \lambda_{8} = \lambda_{16} = \lambda_{17} = \lambda_{18}= \frac{8(2-s_\nu)}{8-s_\nu} \;,
   \\end{eqnarray}
   $$

The non-zero equilibrium moments are defined as

.. math::
   :nowrap:

   $$
   \\begin{eqnarray}
     m_1^{eq} &=& (j_x^2+j_y^2+j_z^2) - \alpha |\textbf{C}|, \\
     m_9^{eq} &=& (2j_x^2-j_y^2-j_z^2)+ \alpha \frac{|\textbf{C}|}{2}(2n_x^2-n_y^2-n_z^2), \\
     m_{11}^{eq} &=& (j_y^2-j_z^2) + \alpha \frac{|\textbf{C}|}{2}(n_y^2-n_z^2), \\
     m_{13}^{eq} &=& j_x j_y + \alpha \frac{|\textbf{C}|}{2} n_x n_y\;, \\
     m_{14}^{eq} &=& j_y j_z + \alpha \frac{|\textbf{C}|}{2} n_y n_z\;, \\
     m_{15}^{eq} &=& j_x j_z + \alpha \frac{|\textbf{C}|}{2} n_x n_z\;, 
   \\end{eqnarray}
   $$


   
****************************
Boundary Conditions
****************************

The following external boundary conditions are supported by ``lbpm_color_simulator``
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

