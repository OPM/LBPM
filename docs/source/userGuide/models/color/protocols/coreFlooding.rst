======================================
Color model -- Core Flooding
======================================

The core flooding protocol is designed to mimic SCAL experiments where one
immiscible fluid is injected into the sample at a constant rate, displacing the
other fluid. The core flooding protocol relies on a flux boundary condition
to ensure that fluid is injected into the sample at a constant rate. The flux
boundary condition implements a time-varying pressure boundary condition that
adapts to ensure a constant volumetric flux. Details for the flux boundary
condition are available
(see: https://doi.org/10.1016/j.compfluid.2020.104670)

.. code-block:: bash

    protocol = "core flooding"
    

To match experimental conditions, it is usually important to match the capillary
number, which is

.. math::
   \mbox{Ca} = \frac{\mu u_z}{\gamma}


where :math:`\mu` is the dynamic viscosity, :math:`u_z` is the fluid
(usually water) velocity and :math:`\gamma` is the interfacial tension.
The volumetric flow rate is related to the fluid velocity based on

.. math::
   Q_z = \epsilon C_{xy} u_z

where :math:`C_{xy}` is the cross-sectional area and :math:`\epsilon`
is the porosity. Given a particular experimental system 
self-similar conditions can be determined for the lattice Boltzmann
simulation by matching the non-dimensional :math:`mbox{Ca}`. It is nearly
awlays advantageous for the timestep to be as large as possible so
that time-to-solution will be more favorable. This is accomplished by

* use a high value for the numerical surface tension (e.g. ``alpha=1.0e-2``)
* use a small value for the fluid viscosity (e.g. ``tau_w = 0.7`` and ``tau_n = 0.7`` )
* determine the volumetric flow rate needed to match :math:`\mbox{Ca}`

For the color LBM the interfacial tension is
:math:`\gamma = 6 \alpha` and the dynamic viscosity is :math:`\mu =  \rho(\tau-1/2)/3`,
where the units are relative to the lattice spacing, timestep and mass
density. Agreemetn between the experimental and simulated values for
:math:`\mbox{Ca}` is ensured by setting the volumetric flux

.. math::
   Q_z = \frac{\epsilon C_{xy} \gamma }{\mu} \mbox{Ca}

where the LB units of the volumetric flux will be voxels per timestep.

In some situations it may also be important to match other non-dimensional numbers,
such as the viscosity ratio, density ratio, and/or Ohnesorge/Laplace number. This
can be accomplished with an analogous procedure. Enforcing additional constraints
will necessarily restrict the LB parameters that can be used, which are ultimately
manifested as a constraint on the size of a timestep. 




