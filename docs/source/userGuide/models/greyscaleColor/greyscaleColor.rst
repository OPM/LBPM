=============================================
Greyscale Color Model
=============================================

The LBPM greyscale lattice Boltzmann model is constructed to approximate the
solution of the Darcy-Brinkman equations in grey regions coupled to a color model implementation
solution in open regions. A simple constitutive form is used to model the relative
permeability in the grey regions,


***************************
Parameters
***************************

The essential model parameters for the color model are

- ``alpha`` -- control the interfacial tension between fluids -- :math:`0 < \alpha < 0.01`
- ``beta`` -- control the width of the interface -- :math:`\beta < 1`
- ``tauA`` -- control the viscosity of fluid A -- :math:`0.7 < \tau_A < 1.5`
- ``tauB`` -- control the viscosity of fluid B -- :math:`0.7 < \tau_B < 1.5`
- ``rhoA`` -- control the viscosity of fluid A -- :math:`0.05 < \rho_A < 1.0`
- ``rhoB`` -- control the viscosity of fluid B -- :math:`0.05 < \rho_B < 1.0`
- ``greyscale_endpoint_A`` -- control the endpoint for greyscale components -- :math:`S_b^{r,a}`
- ``greyscale_endpoint_B`` -- control the endpoint for greyscale components -- :math:`S_b^{r,b}`

  
****************************
Formulation
****************************

The greyscale color model extends the conventional two-fluid color model such that flow through
micro-porous "grey" voxels can also be treated. Extensions to the formulation are made to:
(1) implement momentum transport behavior consistent with the single phase greyscale model
(2) adapt the re-coloring term to allow for transport through microporosity.
Mass transport LBEs are constructed to model the behavior for each fluid. Within the open pore
regions, mass transport LBEs recover the behavior of the standard two-fluid color model.
Within grey voxels the standard recoloring term is disabled so that the mass transport behavior
can be described by a different rule within the microporous regions. The endpoints are
specified based on the saturation of fluid B at which fluid A ceases to percolate, :math:`S_b^{r,a}`
and the saturation of fluid B at which fluid B ceases to percolate, :math:`S_b^{r,b}`.
The endpoints should be independently specified for each class of microporosity that is labeled
in the input image. Between the endpoint values, the effective permeability is paramaterized based on the value

.. math::
   :nowrap:

    $$
    S_{ab} = \frac{S_b - S_b^{r,b}}{S_b^{r,a} - S_b^{r,b}}
    $$


At the endpoints, the effective permeability is provided as an input parameter for
each fluid. When :math:`S_b=S_b^{r,b}` the effective permeability of fluid B is
zero and :math:`K_a=K^{r,a}`. Between the endpoints the Corey model predicts the
effective permeability for fluid A according to


.. math::
   :nowrap:

      $$
      K_a = K^{r,a} (1-S_{ab})^{\lambda^a}
      $$

where :math:`\lambda^a=2` is the Corey exponent. Likewise, 
When :math:`S_b=S_b^{r,a}` the effective permeability for fluid A will be zero,
and :math:`K_b=K^{r,b}` with 

.. math::
   :nowrap:

      $$
      K_b = K^{r,b} S_{ab}^{\lambda^b}
      $$

      
Two LBEs are constructed to model the mass transport,

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

The recoloring step incorporates the standard color model
rule to model anti-diffusion at the interface within the open pores. Within grey regions,
the anti-diffusion term must be disabled, since within grey voxels the length scale for fluid features 
is smaller than the interface width produced from the color model. Within grey voxels the
two fluids are permitted to freely mix between the endpoints. Beyond the endpoints, the recoloring
term is used to drive spontaneous imbibition into the grey voxels


.. math::
   :nowrap:

      $$
      R_c = 
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


Due to the nature of the contribution of the porosity to the pressure term in the Chapman-Enskog
expansion, periodic boundary conditions are recommended for ``lbpm_greyscaleColor_simulator``
and can be set by setting the ``BC`` key values in the ``Domain`` section of the
input file database

- ``BC = 0`` -- fully periodic boundary conditions

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

