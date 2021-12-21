=============================================
Greyscale model
=============================================

The LBPM greyscale lattice Boltzmann model is constructed to approximate the
solution of the Darcy-Brinkman equations in grey regions, coupled to a Navier-Stokes
solution in open regions. To use the greyscale model, the input image should be segmented
into "grey" and open regions. Each particular "grey" label should be assigned both
a porosity and permeability value. 

A typical command to launch the LBPM color simulator is as follows

```
mpirun -np $NUMPROCS lbpm_greyscale_simulator input.db
```

where ``$NUMPROCS`` is the number of MPI processors to be used and ``input.db`` is
the name of the input database that provides the simulation parameters.
Note that the specific syntax to launch MPI tasks may vary depending on your system.
For additional details please refer to your local system documentation.

***************************
Model parameters
***************************

The essential model parameters for the single phase greyscale model are

- ``tau`` -- control the fluid viscosity -- :math:`0.7 < \tau < 1.5`

The kinematic viscosity is given by

***************************
Model formulation
***************************

A D3Q19 LBE is constructed to describe the momentum transport

.. math::
   :nowrap:

      $$
      f_q(\bm{x}_i + \bm{\xi}_q \delta t,t + \delta t) - f_q(\bm{x}_i,t) =
      \sum^{Q-1}_{k=0} M^{-1}_{qk} S_{kk} (m_k^{eq}-m_k)  + \sum^{Q-1}_{k=0} M^{-1}_{qk} (1-\frac{S_{kk}}{2}) \hat{F}_q\;,
      $$


The force is imposed based on the construction developed by Guo et al

.. math::
   :nowrap:

      $$
      F_i = \rho_0 \omega_i \left[\frac{\bm{e}_i \cdot \bm{a}}{c_s^2} +
      \frac{\bm{u} \bm{a}:(\bm{e}_i \bm{e}_i -c_s^2 \mathcal{I})}{\epsilon c_s^4}   \right] ,
      $$


where :math:`c_s^2 = 1/3` is the speed of sound for the D3Q19 lattice.
The acceleration includes contributions due to the external driving force :math:`\bm{g}`
as well as a drag force due to the permeability :math:`K` and flow rate :math:`\bm{u}` with the
porosity :math:`\epsilon` and  viscosity :math:`\nu` determining the net forces acting within
a grey voxel

.. math::
   :nowrap:

      $$
      \bm{a} = - \frac{\epsilon \nu}{K} \bm{u} + \bm{F}_{cp}/\rho_0 + \epsilon \bm{g},
      $$

The flow velocity is defined as

.. math::
   :nowrap:

      $$
      \rho_0 \bm{u} = \sum_i \bm{e}_i f_i + \frac{\delta t}{2} \rho_0 \bm{a}.
      $$

Combining the previous expressions, 

.. math::
   :nowrap:

      $$
      \bm{u} = \frac{\frac{1}{\rho_0}\sum_i \bm{e}_i f_i + \frac{\delta t}{2}\epsilon \bm{g} +
      \frac{\delta t}{2} \frac{\bm{F}_{cp}}{\rho_0}}{1+ \frac{\delta t}{2} \frac{\epsilon \nu}{K}}
      $$


