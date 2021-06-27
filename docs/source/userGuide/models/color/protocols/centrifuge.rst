======================================
Color model -- Centrifuge Protocol
======================================

The centrifuge protocol is designed to mimic SCAL centrifuge experiments that
are used to infer the capillary pressure. The LBPM centrifuge protocol is
constructed as an unsteady simulation with constant pressure boundary conditions
and zero pressure drop across the sample. This will enforce the following key values

* ``BC = 3`` -- constant pressure boundary condition
* ``din = 1.0`` -- inlet pressure value 
* ``dout = 1.0`` -- outlet pressure value

By default LBPM will populate the inlet reservoir with fluid A (usually the non-wetting fluid)
and the outlet reservoir with fluid B (usually water). Flow is induced by setting an external
body force to generate displacement in the ``z`` direction. If the body force is set to
zero, e.g. ``F = 0, 0, 0``, the simulation will produce spontaneous imbibition, with the
balancing point being determined based on zero pressure drop across the sample. Setting
an external body force will shift the capillary pressure. Setting a positive force will
cause fluid A to be forced into the sample. Once steady conditions are achieved,
the pressure of fluid A will be larger than fluid B. Alternatively, if the driving force is
negative then fluid B will be forced into the sample, and the steady-state configuration
will stabilize to a configuration where fluid B has a larger pressure compared to fluid A.
The capillary pressure is thereby inferred based on the body force.

In a conventional SCAL experiment the centrifugal forces are proportional to the density
difference between fluids. While this is still true for LBPM simulation, the body force will
still be effective even if there is no difference in density between the fluids.
This is because a positive body force will favor a larger saturation of fluid A
(positive capillary pressure ) whereas a negative body force will favor a lower
saturation of fluid A (negative capillary pressure). 


To enable the ``centrifuge`` protocol such that the pressure of fluid B is higher than
fluid A, the following keys can be set. Increasing the body force will lead to a larger
capillary pressure

.. code-block:: bash
		
    protocol = "centrifuge"
    F = 0, 0, -1.0e-5
    

