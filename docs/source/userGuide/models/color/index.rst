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


Two LBEs are constructed to model the mass transport, incorporating the anti-diffusion

.. math::
   :nowrap:

   $$
   A_q(\bm{x} + \bm{\xi}_q \delta t, t+\delta t) = w_q N_a \Big[1 + \frac{\bm{u} \cdot \bm{\xi}_q}{c_s^2} 
       + \beta  \frac{N_b}{N_a+N_b} \bm{n} \cdot \bm{\xi}_q\Big] \;
   $$

.. math::
   :nowrap:

   $$
   B_q(\bm{x} + \bm{\xi}_q \delta t, t+\delta t) = 
       w_q N_b \Big[1 + \frac{\bm{u} \cdot \bm{\xi}_q}{c_s^2}
       - \beta  \frac{N_a}{N_a+N_b} \bm{n} \cdot \bm{\xi}_q\Big]\;, 
   $$

The number density for each fluid is obtained from the sum of the mass transport distributions

.. math::
   :nowrap:

   $$
   N_a = \sum_q A_q\;, \quad    N_b = \sum_q B_q\; 
   $$

   
The phase indicator field is then defined as 

.. math::
   :nowrap:

   $$
   \phi = \frac{N_a-N_b}{N_a+N_b}
   $$

The fluid density and kinematic viscosity are determined based on linear interpolation

   
.. math::
   :nowrap:

   $$
    \rho_0 = \frac{(1+\phi) \rho_n}{2}+ \frac{(1-\phi) \rho_w}{2} \;,
   $$

.. math::
   :nowrap:

   $$
    s_\nu = \frac{(1+\phi)}{2\tau_n} +\frac{(1-\phi)}{2\tau_w} \;,
   $$

where

.. math::
   :nowrap:

   $$
    \nu_w = \frac{1}{3}\Big(\tau_w - \frac{1}{2} \Big) \;, \quad
    \nu_n = \frac{1}{3}\Big(\tau_n - \frac{1}{2} \Big) \;.
   $$


These values are then used to model the momentum transport.
The LBE governing momentum transport is defined based on a MRT relaxation process with additional
terms to account for the interfacial stresses

.. math::
   :nowrap:

   $$
      f_q(\bm{x}_i + \bm{\xi}_q \delta t,t + \delta t) - f_q(\bm{x}_i,t) = \sum^{Q-1}_{k=0} M^{-1}_{qk} \lambda_{k} (m_k^{eq}-m_k) + w_q \bm{\xi}_q \cdot \frac{\bm{F}}{c_s^2} \;,
   $$

Where :math:`\bm{F}` is an external body force and :math:`c_s^2 = 1/3` is the speed of sound for the LB model.
The moments are linearly indepdendent:

.. math::
   :nowrap:

   $$
      m_k = \sum_{q=0}^{18} M_{qk} f_q\;.
   $$

   
The relaxation parameters are determined from the relaxation time:

.. math::
   :nowrap:

   $$
     \lambda_1 =  \lambda_2=  \lambda_9 = \lambda_{10}= \lambda_{11}= \lambda_{12}= \lambda_{13}= \lambda_{14}= \lambda_{15} = s_\nu \;,
   $$
   
.. math::
   :nowrap:
      
    $$
     \lambda_{4}= \lambda_{6}= \lambda_{8} = \lambda_{16} = \lambda_{17} = \lambda_{18}= \frac{8(2-s_\nu)}{8-s_\nu} \;,
   $$

The non-zero equilibrium moments are defined as

.. math::
   :nowrap:

   $$
     m_1^{eq} = 19\frac{ j_x^2+j_y^2+j_z^2}{\rho_0} - 11\rho - 19 \alpha |\textbf{C}|, \\
   $$     

.. math::
   :nowrap:

   $$
     m_2^{eq} = 3\rho - \frac{11( j_x^2+j_y^2+j_z^2)}{2\rho_0}, \\
   $$     

.. math::
   :nowrap:

   $$
     m_4^{eq} = -\frac{2 j_x}{3}, \\
   $$     

.. math::
   :nowrap:

   $$
     m_6^{eq} = -\frac{2 j_y}{3}, \\
   $$     

.. math::
   :nowrap:

   $$
     m_8^{eq} = -\frac{2 j_z}{3}, \\
   $$     

.. math::
   :nowrap:

   $$     
     m_9^{eq} = \frac{2j_x^2-j_y^2-j_z^2}{\rho_0}+ \alpha \frac{|\textbf{C}|}{2}(2n_x^2-n_y^2-n_z^2), \\
   $$     

.. math::
   :nowrap:

   $$     
     m_{11}^{eq} = \frac{j_y^2-j_z^2}{\rho_0} + \alpha \frac{|\textbf{C}|}{2}(n_y^2-n_z^2), \\
   $$     

.. math::
   :nowrap:

   $$     
     m_{13}^{eq} = \frac{j_x j_y}{\rho_0} + \alpha \frac{|\textbf{C}|}{2} n_x n_y\;, \\
   $$     

.. math::
   :nowrap:

   $$     
     m_{14}^{eq} = \frac{j_y j_z}{\rho_0} + \alpha \frac{|\textbf{C}|}{2} n_y n_z\;, \\
   $$     

.. math::
   :nowrap:

   $$     
     m_{15}^{eq} = \frac{j_x j_z}{\rho_0} + \alpha \frac{|\textbf{C}|}{2} n_x n_z\;. 
   $$

where the color gradient is determined from the phase indicator field

.. math::
   :nowrap:

   $$
   \textbf{C}=\nabla \phi\;.
   $$

and the unit normal vector is

.. math::
   :nowrap:

   $$
     \bm{n} = \frac{\textbf{C}}{|\textbf{C}|}\;.
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

****************
Example data
****************

Example data can be downloaded from https://www.digitalrocksportal.org/projects/326 
