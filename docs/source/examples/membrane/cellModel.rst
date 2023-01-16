********************************
Membrane Charging Dynamics
********************************

In this example, we consider membrane charging dynamics for a simple cell. 

For the case considered in ``example/SingleCell`` an input membrane geometry is provided in the
file ``Bacterium.swc``, which specifies an oblong cell shape, relying on the ``.swc`` file format that
is commonly used to approximate neuron structures. The case considered is the four ion membrane transport
problem considered in Figure 4 from McClure & Li 

The cell simulation is performed by the executable ``lbpm_nernst_planck_cell_simulator``, which is launched
in the same way as other parallel tools

.. code:: bash

	  mpirun -np 2 $LBPM_BIN/lbpm_nernst_planck_cell_simulator Bacterium.db


The input file ``Bacterium.db`` specifies the following

.. code:: c

	  MultiphysController {
	     timestepMax = 25000
	     num_iter_Ion_List = 4
	     analysis_interval  = 100
	     tolerance = 1.0e-9
	     visualization_interval = 1000        // Frequency to write visualization data
	  }
	  Ions {
	     use_membrane = true
	     Restart = false
	     MembraneIonConcentrationList = 150.0e-3, 10.0e-3, 15.0e-3, 155.0e-3 //user-input unit: [mol/m^3]
	     temperature = 293.15 //unit [K]
	     number_ion_species = 4  //number of ions
	     tauList = 1.0, 1.0, 1.0, 1.0
	     IonDiffusivityList = 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9 //user-input unit: [m^2/sec]
	     IonValenceList = 1, -1, 1, -1 //valence charge of ions; dimensionless; positive/negative integer
	     IonConcentrationList = 4.0e-3, 20.0e-3, 16.0e-3, 0.0e-3 //user-input unit: [mol/m^3]
	     BC_Solid = 0 //solid boundary condition; 0=non-flux BC; 1=surface ion concentration
	     FluidVelDummy = 0.0, 0.0, 0.0 // dummy fluid velocity for debugging
	     BC_InletList = 0, 0, 0, 0
	     BC_OutletList = 0, 0, 0, 0
	  }
	  Poisson {
  	     lattice_scheme = "D3Q19"
	     epsilonR = 78.5 //fluid dielectric constant [dimensionless]
	     BC_Inlet  = 0  // ->1: fixed electric potential; ->2: sine/cosine periodic electric potential
	     BC_Outlet = 0  // ->1: fixed electric potential; ->2: sine/cosine periodic electric potential
	     //--------------------------------------------------------------------------
	     //--------------------------------------------------------------------------
	     BC_Solid = 2 //solid boundary condition; 1=surface potential; 2=surface charge density
	     SolidLabels = 0 //solid labels for assigning solid boundary condition
	     SolidValues = 0 //if surface potential, unit=[V]; if surface charge density, unit=[C/m^2]
	     WriteLog = true //write convergence log for LB-Poisson solver
	     //------------------------------ advanced setting ------------------------------------
	     timestepMax = 4000 //max timestep for obtaining steady-state electrical potential
	     analysis_interval  = 25 //timestep checking steady-state convergence
	     tolerance = 1.0e-10  //stopping criterion for steady-state solution
	     InitialValueLabels = 1, 2
	     InitialValues = 0.0, 0.0
	  }
	  Domain {
  	     Filename = "Bacterium.swc"
	     nproc = 2, 1, 1     // Number of processors (Npx,Npy,Npz)
	     n = 64, 64, 64      // Size of local domain (Nx,Ny,Nz)
	     N = 128, 64, 64         // size of the input image
	     voxel_length = 0.01   //resolution; user-input unit: [um]
	     BC = 0              // Boundary condition type
	     ReadType = "swc"
	     ReadValues  = 0, 1, 2
	     WriteValues = 0, 1, 2
	  }
	  Analysis {
	     analysis_interval = 100
	     subphase_analysis_interval = 50    // Frequency to perform analysis
	     restart_interval = 5000    // Frequency to write restart data
	     restart_file = "Restart"    // Filename to use for restart file (will append rank)
	     N_threads    = 4            // Number of threads to use
	     load_balance = "independent" // Load balance method to use: "none", "default", "independent"
	  }
	  Visualization {
	     save_electric_potential = true
	     save_concentration = true
	     save_velocity = false
	  }
	  Membrane {
	     MembraneLabels = 2
	     VoltageThreshold = 0.0, 0.0, 0.0, 0.0
	     MassFractionIn = 1e-1, 1.0, 5e-3, 0.0
	     MassFractionOut = 1e-1, 1.0, 5e-3, 0.0
	     ThresholdMassFractionIn = 1e-1, 1.0, 5e-3, 0.0
	     ThresholdMassFractionOut = 1e-1, 1.0, 5e-3, 0.0
	  }
	  
*******************
Example Output
*******************

Successful output looks like the following

.. code:: bash

	  ********************************************************
	  Running LBPM Nernst-Planck Membrane solver 
	  ********************************************************
	  .... Read membrane permeability (MassFractionIn) 
	  .... Read membrane permeability (MassFractionOut) 
	  .... Read membrane permeability (ThresholdMassFractionIn) 
	  .... Read membrane permeability (ThresholdMassFractionOut) 
	  .... Read MembraneIonConcentrationList 
	  voxel length = 0.010000 micron 
	  voxel length = 0.010000 micron 
	  Reading SWC file...
	      Number of lines in SWC file: 7
	     Number of lines extracted is: 7
	     shift swc data by 0.150000, 0.140000, 0.140000 
	  Media porosity = 1.000000 
	  LB Ion Solver: Initialized solid phase & converting to Signed Distance function 
	      Domain set.
	  LB Ion Solver: Create ScaLBL_Communicator 
	  LB Ion Solver: Set up memory efficient layout 
	  LB Ion Solver: Allocating distributions 
	  LB Ion Solver: Setting up device map and neighbor list 
	  **** Creating membrane data structure ****** 
	     Number of active lattice sites (rank = 0): 262160 
	     Membrane labels: 1 
		label=2, volume fraction = 0.133917
	  Creating membrane data structure...
	     Copy initial neighborlist... 
	     Cut membrane links... 
	     (cut 7105 links crossing membrane) 
	     Construct membrane data structures... 
	     Create device data structures... 
	     Construct communication data structures... 
	  Ion model setup complete
	  Analyze system with sub-domain size = 66 x 66 x 66 
	  Set up analysis routines for 4 ions 
	  LB Ion Solver: initializing D3Q7 distributions
	     ...initializing based on membrane list 
	  .... Set concentration(0): inside=0.15 [mol/m^3], outside=0.004 [mol/m^3] 
	  .... Set concentration(1): inside=0.01 [mol/m^3], outside=0.02 [mol/m^3] 
	  .... Set concentration(2): inside=0.015 [mol/m^3], outside=0.016 [mol/m^3] 
	  .... Set concentration(3): inside=0.155 [mol/m^3], outside=0 [mol/m^3] 
	  LB Ion Solver: initializing charge density
	  LB Ion Solver: solid boundary: non-flux boundary is assigned
	  LB Ion Solver: inlet boundary for Ion 1 is periodic 
	  LB Ion Solver: outlet boundary for Ion 1 is periodic 
	  LB Ion Solver: inlet boundary for Ion 2 is periodic 
	  LB Ion Solver: outlet boundary for Ion 2 is periodic 
	  LB Ion Solver: inlet boundary for Ion 3 is periodic 
	  LB Ion Solver: outlet boundary for Ion 3 is periodic 
	  LB Ion Solver: inlet boundary for Ion 4 is periodic 
	  LB Ion Solver: outlet boundary for Ion 4 is periodic 
	  *****************************************************
	  LB Ion Transport Solver: 
		Ion 1: LB relaxation tau = 1
			Time conversion factor: 1.25e-08 [sec/lt]
			Internal iteration: 2 [lt]
		Ion 2: LB relaxation tau = 1
			Time conversion factor: 1.25e-08 [sec/lt]
			Internal iteration: 2 [lt]
		Ion 3: LB relaxation tau = 1
			Time conversion factor: 1.25e-08 [sec/lt]
			Internal iteration: 2 [lt]
		Ion 4: LB relaxation tau = 1
			Time conversion factor: 1.25e-08 [sec/lt]
			Internal iteration: 2 [lt]
	  *****************************************************
	  Ion model initialized 
	  Main loop time_conv computed from ion 1: 2.5e-08[s/lt]
	  Main loop time_conv computed from ion 2: 2.5e-08[s/lt]
	  Main loop time_conv computed from ion 3: 2.5e-08[s/lt]
	  Main loop time_conv computed from ion 4: 2.5e-08[s/lt]
	  ***********************************************************************************
	  LB-Poisson Solver: steady-state MaxTimeStep = 4000; steady-state tolerance = 1e-10 
			     LB relaxation tau = 3.5 
	  ***********************************************************************************
	  LB-Poisson Solver: Use averaged MSE to check solution convergence.
	  LB-Poisson Solver: Use D3Q19 lattice structure.
	  voxel length = 0.010000 micron 
	  voxel length = 0.010000 micron 
	  Reading SWC file...
	      Number of lines in SWC file: 7
	     Number of lines extracted is: 7
	     shift swc data by 0.150000, 0.140000, 0.140000 
	  Media porosity = 1.000000 
	  LB-Poisson Solver: Initialized solid phase & converting to Signed Distance function 
	      Domain set.
	  LB-Poisson Solver: Create ScaLBL_Communicator 
	  LB-Poisson Solver: Set up memory efficient layout 
	  LB-Poisson Solver: Allocating distributions 
	  LB-Poisson Solver: Setting up device map and neighbor list 
	   .... LB-Poisson Solver: check  neighbor list 
	   .... LB-Poisson Solver: copy  neighbor list to GPU 
	  Poisson solver created 
	  LB-Poisson Solver: initializing D3Q19 distributions
	  LB-Poisson Solver: number of Poisson solid labels: 1 
	     label=0, surface potential=0 [V], volume fraction=0
	  LB-Poisson Solver: number of Poisson initial-value labels: 2 
	     label=1, initial potential=0 [V], volume fraction=0.96
	     label=2, initial potential=0 [V], volume fraction=0.13
	     POISSON MODEL: Reading restart file! 
	  Poisson solver initialized 
	     ... getting Poisson solver error 
	  -------------------------------------------------------------------
	  set coefficients 
	  ********************************************************
	  CPU time = 0.008526 
	  Lattice update rate (per core)= 30.749833 MLUPS 
	  Lattice update rate (total)= 61.499666 MLUPS 
	  ********************************************************
