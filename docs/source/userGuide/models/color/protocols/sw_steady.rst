======================================
Color model -- Sw Steady
======================================

The water saturation steady state protocol is identical to the centrifuge protocol, with
the exception that the simulation explicity converges on the saturation state of fluid
A.

That is, the simulation exits when

.. math::
   :nowrap:

   $$
   \frac{\left | S_{w, i+1} - S_{w, i} \right |}{S_{w, i}} \le \epsilon
   $$

averaged over the ``analysis_interval`` where :math:`S_{w,i}` is the saturation of fluid
A at step :math:`i` and :math:`\epsilon` is an allowed convergence threshold, or when
the ``timestepMax`` is reached, whichever occurs first.

By default, :math:`\epsilon` is set to zero (i.e., the simulation continues through the
``timestepMax``), but can be set via ``tolerance`` within the ``Analysis`` section of
the input file database as shown below.


.. code-block:: c

   Color {
      protocol = "sw_steady"
      timestepMax = 1000000              // maximum timtestep
      alpha = 0.005                      // controls interfacial tension
      rhoA = 1.0                         // controls the density of fluid A
      rhoB = 1.0                         // controls the density of fluid B
      tauA = 0.7                         // controls the viscosity of fluid A
      tauB = 0.7                         // controls the viscosity of fluid B 
      F = 0, 0, -1.0e-5                  // body force
      din = 1.0                          // inlet density (controls pressure)
      dout = 1.0                         // outlet density (controls pressure)   
      WettingConvention = "SCAL"         // convention for sign of wetting affinity
      ComponentLabels = 0, -1, -2        // image labels for solid voxels
      ComponentAffinity = 1.0, 1.0, 0.6  // controls the wetting affinity for each label
      Restart = false
   }
   Domain {
      Filename = "Bentheimer_LB_sim_intermediate_oil_wet_Sw_0p37.raw"  
      ReadType = "8bit"              // data type
      N = 900, 900, 1600             // size of original image
      nproc = 2, 2, 2                // process grid
      n = 200, 200, 200              // sub-domain size
      offset = 300, 300, 300         // offset to read sub-domain
      voxel_length = 1.66            // voxel length (in microns)
      ReadValues = -2, -1, 0, 1, 2   // labels within the original image
      WriteValues = -2, -1, 0, 1, 2  // associated labels to be used by LBPM
      BC = 3                         // boundary condition type (0 for periodic)
   }
   Analysis {
      analysis_interval = 1000           // logging interval for timelog.csv
      subphase_analysis_interval = 5000  // loggging interval for subphase.csv
      visualization_interval = 100000    // interval to write visualization files
      N_threads = 4                      // number of analysis threads (GPU version only)
      restart_interval = 1000000         // interval to write restart file
      restart_file = "Restart"           // base name of restart file
      tolerance = 1e-9                   // Sw convergence tolerance
   }
   Visualization {
      write_silo = true     // write SILO databases with assigned variables
      save_8bit_raw = true  // write labeled 8-bit binary files with phase assignments
      save_phase_field = true  // save phase field within SILO database
      save_pressure = false    // save pressure field within SILO database
      save_velocity = false    // save velocity field within SILO database
   }
