********************************
Steady-state flow (color model)
********************************

In this example we simulate a steady-state flow with a constant driving force. This will enforce a periodic boundary condition
in all directions. While the driving force may be set in any direction, we will set it in the z-direction to be consistent
with the convention for pressure and velocity boundary conditions. 


For the case considered in ``example/DiscPack`` we specify the following information in the input file

.. code:: c

	Domain {
	   Filename = "discs_3x128x128.raw.morphdrain.raw"
	   ReadType = "8bit"    // data type
	   N = 3, 128, 128       // size of original image
	   nproc = 1, 2, 2       // process grid
	   n = 3, 64, 64         // sub-domain size
	   voxel_length = 1.0    // voxel length (in microns)
	   ReadValues = 0, 1, 2  // labels within the original image
	   WriteValues = 0, 1, 2 // associated labels to be used by LBPM
	   BC = 0                // fully periodic BC
	   Sw = 0.35             // target saturation for morphological tools
	}

	Color {
	   protocol = "fractional flow"
	   capillary_number = 1e-5            // capillary number for the displacement, positive="oil injection"
	   timestepMax = 500000               // maximum timtestep
	   alpha = 0.005                      // controls interfacial tension
	   rhoA = 1.0                         // controls the density of fluid A
	   rhoB = 1.0                         // controls the density of fluid B
	   tauA = 0.7                         // controls the viscosity of fluid A
	   tauB = 0.7                         // controls the viscosity of fluid B
	   F = 0, 0, 1e-5                     // body force
	   WettingConvention = "SCAL"
	   ComponentLabels = 0        // image labels for solid voxels
	   ComponentAffinity = 0.9  // controls the wetting affinity for each label
	   Restart = false
	}
	Analysis {
	   analysis_interval = 1000           // logging interval for timelog.csv
	   subphase_analysis_interval = 500000  // loggging interval for subphase.csv
	   N_threads = 4                      // number of analysis threads (GPU version only)
	   visualization_interval = 1000000    // interval to write visualization files
	   restart_interval = 10000000         // interval to write restart file
	   restart_file = "Restart"           // base name of restart file
	}
	Visualization {
	   format = "hdf5"
	   write_silo = true        // write SILO databases with assigned variables
	   save_8bit_raw = true      // write labeled 8-bit binary files with phase assignments
	   save_phase_field = true  // save phase field within SILO database
	   save_pressure = true     // save pressure field within SILO database
	   save_velocity = false     // save velocity field within SILO database
	}
	FlowAdaptor {
	   min_steady_timesteps = 250000       // minimum number of timesteps per steady point
	   max_steady_timesteps = 300000       // maximum number of timesteps per steady point
	   fractional_flow_increment = 0.1   // parameter that controls rate of mass seeding
	   skip_timesteps = 10000             // number of timesteps to spend in flow adaptor
	   endpoint_threshold = 0.1           // endpoint exit criterion
	}

	  
Once this has been set, we launch ``lbpm_color_simulator`` in the same way as other parallel tools

.. code:: bash

	  mpirun -np 4 $LBPM_BIN/lbpm_color_simulator input.db

Successful output looks like the following


.. code:: bash

      ********************************************************
      Running Color LBM	
      ********************************************************
      voxel length = 1.000000 micron 
      voxel length = 1.000000 micron 
      Input media: discs_3x128x128.raw.morphdrain.raw
      Relabeling 3 values
      oldvalue=0, newvalue =0 
      oldvalue=1, newvalue =1 
      oldvalue=2, newvalue =2 
      Dimensions of segmented image: 3 x 128 x 128 
      Reading 8-bit input data 
      Read segmented data from discs_3x128x128.raw.morphdrain.raw 
      Label=0, Count=11862 
      Label=1, Count=26430 
      Label=2, Count=10860 
      Distributing subdomains across 4 processors 
      Process grid: 1 x 2 x 2 
      Subdomain size: 3 x 64 x 64 
      Size of transition region: 0 
      Media porosity = 0.758667 
      Initialized solid phase -- Converting to Signed Distance function 
      Domain set.
      Create ScaLBL_Communicator 
      Set up memory efficient layout, 9090 | 9120 | 21780 
      Allocating distributions 
      Setting up device map and neighbor list 
      Component labels: 1 
	 label=0, affinity=-0.900000, volume fraction==0.417582
      Initializing distributions 
      Initializing phase field 
      Affinities - rank 0:
      Main: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 1: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 2: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 3: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 4: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Affinities - rank 0:
      Main: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 1: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 2: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 3: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      Thread 4: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
      ********************************************************
      CPU time = 0.001501 
      Lattice update rate (per core)= 6.074861 MLUPS 
      Lattice update rate (per MPI process)= 6.074861 MLUPS 
	 (flatten density field)  

