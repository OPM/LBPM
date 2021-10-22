======================================
Color model -- Image Sequence Protocol
======================================

The image sequence protocol is designed to perform a set steady-state
simulations based on a sequence of 3D (8-bit) images provided by the user.
The images might be the output of a previous LBPM simulation, a sequence of
(segmented) experimental data, or data generated from a custom routine.
The image sequence protocol will apply the same set of flow conditions
to all images in the sequence. This means

* the image labels and any associated properties are the same
* the external boundary conditions are the same
* the physical simulation parameters are the same

The image sequence protocol does not set boundary conditions by default.
It is up to the user to determine the flow condition, with the understanding
that the same set of will be applied to each image in the sequence.

To enable the image sequence protocol, the following keys should be set
within the ``Color`` section of the input database

.. code-block:: bash
  
   Color {
      protocol = "image sequence"
      image_sequence = "Bentheimer_LB_RelPerm_intermediate_oil_wet_Sw_0p37.raw", "Bentheimer_LB_RelPerm_intermediate_oil_wet_Sw_0p72.raw"
                                         // comma separated list of image names
      capillary_number = 1e-4            // capillary number for the displacement
      timestepMax = 1000000              // maximum timtestep
      alpha = 0.005                      // controls interfacial tension
      rhoA = 1.0                         // controls the density of fluid A
      rhoB = 1.0                         // controls the density of fluid B
      tauA = 0.7                         // controls the viscosity of fluid A
      tauB = 0.7                         // controls the viscosity of fluid B 
      F = 0, 0, 0                        // body force
      WettingConvention = "SCAL"         // convention for sign of wetting affinity
      ComponentLabels = 0, -1, -2        // image labels for solid voxels
      ComponentAffinity = 1.0, 1.0, 0.6  // controls the wetting affinity for each label
      Restart = false
   }
   Domain {
      Filename = "Bentheimer_LB_RelPerm_intermediate_oil_wet_Sw_0p37.raw"  
      ReadType = "8bit"              // data type
      N = 900, 900, 1600             // size of original image
      nproc = 2, 2, 2                // process grid
      n = 200, 200, 200              // sub-domain size
      offset = 300, 300, 300         // offset to read sub-domain
      voxel_length = 1.66            // voxel length (in microns)
      ReadValues = -2, -1, 0, 1, 2   // labels within the original image
      WriteValues = -2, -1, 0, 1, 2  // associated labels to be used by LBPM
      BC = 0                         // boundary condition type (0 for periodic)
   }
   Analysis {
      analysis_interval = 1000           // logging interval for timelog.csv
      subphase_analysis_interval = 5000  // loggging interval for subphase.csv
      visualization_interval = 100000    // interval to write visualization files
      N_threads = 4                      // number of analysis threads (GPU version only)
      restart_interval = 1000000         // interval to write restart file
      restart_file = "Restart"           // base name of restart file
   }
   Visualization {
      write_silo = true     // write SILO databases with assigned variables
      save_8bit_raw = true  // write labeled 8-bit binary files with phase assignments
      save_phase_field = true  // save phase field within SILO database
      save_pressure = false    // save pressure field within SILO database
      save_velocity = false    // save velocity field within SILO database
   }
   FlowAdaptor {
      min_steady_timesteps = 100000     // minimum number of timesteps per steady point
      max_steady_timesteps = 250000     // maximum number of timesteps per steady point
   }

