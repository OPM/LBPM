=============================================
Gauss's Law Model
=============================================

The LBPM Gauss's law solver is designed to solve for the electric field in an ionic fluid. 

.. math::
   :nowrap:

   $$
    \nabla^2_{fe} \psi (\mathbf{x}_i) = \frac{1}{6 \Delta x^2}
    \Bigg( 2 \sum_{q=1}^{6} \psi(\mathbf{x}_i + \bm{\xi}_q \Delta t) 
      +  \sum_{q=7}^{18} \psi(\mathbf{x}_i + \bm{\xi}_q \Delta t)
   - 24 \psi (\mathbf{x}_i) \Bigg) \;,
    $$

The equilibrium functions are defined as

.. math::
   :nowrap:

   $$
    g_q^{eq} =  w_q \psi\;,
   $$

where :math:`w_0=1/2`, :math:`w_q=1/24` for :math:`q=1,\ldots,6` and :math:`w_q=1/48` for :math:`q=7,\ldots,18`

which implies that 

.. math::
   :nowrap:

   $$
    \psi = \sum_{q=0}^{Q} g_q^{eq}\;.
    $$
    
Given a particular initial condition for :math:`\psi`, let us consider application of the standard D3Q19 streaming step based on the equilibrium distributions

.. math::
   :nowrap:

   $$
    g_q^\prime(\mathbf{x}, t) = g_q^{eq}(\mathbf{x}-\bm{\xi}_q\Delta t, t+ \Delta t)\;.
   $$
   
Relative to the solution of Gauss's law, the error is given by

.. math::
   :nowrap:

   $$
   \varepsilon_{\psi} = 
   8 \Big[ -g_0 +  \sum_{q=1}^Q g_q^\prime(\mathbf{x}, t) \Big] 
   + \frac{\rho_e}{\epsilon_r \epsilon_0} \;.
   $$
     
Using the fact that :math:`f_0 = W_0 \psi`, we can compute the value 
:math:`\psi^\prime` that would kill the error. We set :math:`\varepsilon_{\psi}=0`
and rearrange terms to obtain

.. math::
   :nowrap:

   $$
   \psi^\prime (\mathbf{x},t) = \frac{1}{W_0}\Big[   \sum_{q=1}^Q g_q^\prime(\mathbf{x}, t) 
   + \frac{1}{8}\frac{\rho_e}{\epsilon_r \epsilon_0}\Big]  \;.
   $$

The local value of the potential is then updated based on a relaxation scheme, which is controlled by the relaxation time :math:`\tau_\psi`

.. math::
   :nowrap:

   $$
   \psi(\mathbf{x},t+\Delta t) \leftarrow \Big(1 - \frac{1}{\tau_\psi} \Big )\psi (\mathbf{x},t)
   + \frac{1}{\tau_\psi} \psi^\prime (\mathbf{x},t)\;.
   $$
   
The algorithm can then proceed to the next timestep.

****************************
Example Input File
****************************

.. code-block:: c

     MultiphysController {
	 timestepMax = 25000
	 num_iter_Ion_List = 4
	 analysis_interval  = 100
	 tolerance = 1.0e-9
	 visualization_interval = 1000        // Frequency to write visualization data
     }
     Stokes {
	 tau = 1.0
	 F = 0, 0, 0
	 ElectricField = 0, 0, 0 //body electric field; user-input unit: [V/m]
	 nu_phys = 0.889e-6      //fluid kinematic viscosity; user-input unit: [m^2/sec]
     }
     Ions {
	 MembraneIonConcentrationList = 150.0e-3, 10.0e-3, 15.0e-3, 155.0e-3 //user-input unit: [mol/m^3]
	 temperature = 293.15 //unit [K]
	 number_ion_species = 4  //number of ions
	 tauList = 1.0, 1.0, 1.0, 1.0
	 IonDiffusivityList = 1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9 //user-input unit: [m^2/sec]
	 IonValenceList = 1, -1, 1, -1 //valence charge of ions; dimensionless; positive/negative integer
	 IonConcentrationList = 4.0e-3, 20.0e-3, 16.0e-3, 0.0e-3 //user-input unit: [mol/m^3]
	 BC_Solid = 0 //solid boundary condition; 0=non-flux BC; 1=surface ion concentration
	 //SolidLabels = 0 //solid labels for assigning solid boundary condition; ONLY for BC_Solid=1
	 //SolidValues = 1.0e-5 // user-input surface ion concentration unit: [mol/m^2]; ONLY for BC_Solid=1
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
	 // ------------------------------- Testing Utilities ----------------------------------------
	 // ONLY for code debugging; the followings test sine/cosine voltage BCs; disabled by default
	 TestPeriodic = false
	 TestPeriodicTime = 1.0 //unit:[sec]
	 TestPeriodicTimeConv = 0.01 //unit:[sec]
	 TestPeriodicSaveInterval = 0.2 //unit:[sec]
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
