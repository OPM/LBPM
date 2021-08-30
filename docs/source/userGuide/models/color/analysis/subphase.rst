======================================
Color model -- Subphase Analysis
======================================

The subphase analysis routine for the LBPM color model logs a time series
of averaged information to the space-delimited CSV file ``subphase.csv``.
The ``subphase_analysis_interval`` key should be specified to control the interval at
to perform basic analysis, which corresponds to the number of simulation timesteps
to perform between analyses. The subphase analys routine performs the following functions

* analyzes the connectivity of fluid phases using a connected components algorithm
* constructs iso-surfaces to represent interfaces within the system
* computes averages of physical quantities based on the defined entities

Since it is more computationally expensive to carry out these operations compared to the
basic analysis, it may be useful to choose ``subphase_analysis_interval`` to be larger than
``analysis_interval``. Nevertheless, since analysis is performed in memory, it is orders of
magnitude faster than to analyze data *in situ* rather than writing fields to disc. Monitoring
the performance in MLUPS provides an easy way to tune the analysis interval so that the
overall simulation performance is not subject to significant penalties. 

There are three main entities that are constructed for subphase analysis

* :math:`\Omega_w:\phi<0` : region of space filled by fluid w
* :math:`\Omega_n:\phi<0` : region of space filled by fluid n
* :math:`\Omega_i: |\nabla \phi|<0 > \epsilon` : region of space for the interface region

The phase regions are defined identical to what is reported in ``timelog.csv``.
The interface region is defined explicitly as the portion of space where
significant composition gradients are present. For each entity, sub-entities are
constructed by running a connected components algorithm on the set. This is done to
separate the connected part of the phase from the disconnected part. This subset operation
is performed by identifying the largest connected component (based on volume)
and denoting this as the connected part of that region. The remaining portion of the
phase is lumped into a disconnected sub-entity. Subphase analysis is therefore performed
for the following six entities

* ``wc`` -- connected part of phase w
* ``wd`` -- disconnected part of phase w
* ``nc`` -- connected part of phase n
* ``nd`` -- disconnected part of phase n
* ``ic`` -- connected part of interface region
* ``id`` -- disconnected part of interface region

For each entity :math:`\Omega_k` with :math:`k\in\{wc,wd,nc,nd,ic,id\}`
an isosurface is constructed to approximate the boundary of the region,
:math:`\Gamma_k`. Once each region is identified, the following measures are determined

**Geometric invariants**

* Volume -- :math:`V_k=\int_{\Omega_k} dV`
* Surface area -- :math:`A_k=\int_{\Gamma_k} dS`
* Integral mean curvature -- :math:`H_k=\int_{\Gamma_k} \frac 12 (\kappa_1 + \kappa_2) dS`
* Euler characteristic -- :math:`\chi_k= \frac{1}{4\pi} \int_{\Gamma_k} \kappa_1 \kappa_2 dS`
      
**Conserved quantities**

* Total mass within the region :math:`M_k=\int_{\Omega_k} \rho dV`
* Total momentum within the region :math:`\mathbf{P}_k=\int_{\Omega_k} \rho \mathbf{u} dV`
* Total kinetic energy within the region :math:`K_k=\frac 12 \int_{\Omega_k} \rho \mathbf{u} \cdot \mathbf{u} dV`

**Thermodynamic quantities**
      
* Pressure -- :math:`p_k=\frac{1}{V_k}\int_{\Omega_k} p dV`
* Solid wetting energy -- :math:`\gamma_s=\int_{\Gamma_s}\gamma dS`
* Viscous dissipation -- :math:`\Phi_k=\int_{\Omega_k}\bm{\varepsilon} : \nabla \mathbf{u} dV`

The total solid wetting energy is determined by integrating the interfacial stresses in the
immediate vicinity of the solid surface :math:`\Gamma_s`. The integral of the
dissipation function is determined based on the viscous stress tensor, denoted by :math:`\bm{\varepsilon}`.


The full list of measures are associated with the labels in ``subphase.csv``

* ``time`` -- timestep 
* ``rn`` -- density for phase n (input parameter)
* ``rw`` -- density for phase w (input parameter)
* ``nun`` -- kinematic viscosity for phase n (input parameter)
* ``nuw`` -- kinematic viscosity for phase w (input parameter)
* ``Fx`` -- external force in x direction (input parameter)
* ``Fy`` -- external force in y direction (input parameter)
* ``Fz`` -- external force in z direction (input parameter)
* ``iftwn`` -- interfacial tension (input parameter)
* ``wet`` -- total solid wetting energy 
* ``pwc`` -- average pressure for connected part of phase w
* ``pwd`` -- average pressure for disconnected part of phase w
* ``pnc`` -- average pressure for connected part of phase n
* ``pnd`` -- average pressure for disconnected part of phase n
* ``Mwc`` -- mass for connected part of phase w 
* ``Mwd`` --mass for disconnected part of phase w
* ``Mwi`` -- mass for phase within diffuse interface region 
* ``Mnc`` -- mass for connected part of phase n 
* ``Mnd`` -- mass for disconnected part of phase n 
* ``Mni`` -- mass for phase n within diffuse interface region
* ``Msw`` -- mass for component w within 2 voxels of solid
* ``Msn`` -- mass for component n within 2 voxels of solid
* ``Pwc_x`` -- x- momentum for connected part of phase w
* ``Pwd_x`` -- x- momentum for disconnected part of phase w
* ``Pwi_x`` -- x- momentum for phase w within diffuse interface
* ``Pnc_x`` -- x- momentum for connected part of phase n
* ``Pnd_x`` -- x- momentum for disconnected part of phase n
* ``Pni_x`` -- x- momentum for phase n within diffuse interface
* ``Psw_x`` -- x- momentum for phase w within 2 voxels of solid
* ``Psn_x`` -- x- momentum for phase n within 2 voxels of solid
* ``Pwc_y`` -- y- momentum for connected part of phase w
* ``Pwd_y`` -- y- momentum for disconnected part of phase w
* ``Pwi_y`` -- y- momentum for phase w within diffuse interface
* ``Pnc_y`` -- y- momentum for connected part of phase n
* ``Pnd_y`` -- y- momentum for connected part of phase n
* ``Pni_y`` -- y- momentum for phase n within diffuse interface
* ``Psw_y`` -- y- momentum for phase w within 2 voxels of solid
* ``Psn_y`` -- y- momentum for phase n within 2 voxels of solid
* ``Pwc_z`` -- z- momentum for connected part of phase w
* ``Pwd_z`` -- z- momentum for disconnected part of phase w
* ``Pwi_z`` -- z- momentum for phase w within diffuse interface
* ``Pnc_z`` -- z- momentum for connected part of phase n
* ``Pnd_z`` -- z- momentum for disconnected part of phase n
* ``Pni_z`` -- z- momentum for phase n within diffuse interface
* ``Psw_z`` -- z- momentum for phase w within 2 voxels of solid
* ``Psn_z`` -- z- momentum for phase n within 2 voxels of solid
* ``Kwc`` -- Kinetic energy for transport within connected part of phase w
* ``Kwd`` -- Kinetic energy for transport within disconnected part of phase w
* ``Kwi`` -- Kinetic energy for transport of phase w within diffuse interface region
* ``Knc`` -- Kinetic energy for transport in connected part of phase n
* ``Knd`` -- Kinetic energy for transport within disconnected part of phase n
* ``Kni`` -- Kinetic energy for transport of phase n within diffuse interface region
* ``Dwc`` -- Viscous dissipation for conneced pathway for phase w
* ``Dwd`` -- Viscous dissipation for disconnected part of phase w
* ``Dnc`` -- Viscous dissipation for connected pathway for phase n
* ``Dnd`` -- Viscous dissipation for disconnected part of phase n
* ``Vwc`` -- Volume for connected pathway for phase w
* ``Awc`` -- Surface area for connected pathway for phase w
* ``Hwc`` -- Integral mean curvature for connected pathway for phase w
* ``Xwc`` -- Euler characteristic for connected pathway for phase w
* ``Vwd`` -- Volume for disconnected phase w
* ``Awd`` -- Surface area for disconnected phase w
* ``Hwd`` -- Integral mean curvature for disconnected phase w
* ``Xwd`` -- Euler characteristic for disconnected phase w
* ``Nwd`` -- Number of connected components in disconnected phase w
* ``Vnc`` -- Volume for connected pathway for phase n
* ``Anc`` -- Surface area for connected pathway for phase n
* ``Hnc`` -- Integral mean curvature for connected pathway for phase n
* ``Xnc`` -- Euler characteristic for connected pathway for phase n
* ``Vnd`` -- Volume for disconnected phase n
* ``And`` -- Surface area for disconnected phase n
* ``Hnd`` -- Integral mean curvature for disconnected phase n
* ``Xnd`` -- Euler characteristic for disconnected phase n
* ``Nnd`` -- number of connected components within disconnected phase n
* ``Vi`` -- volume for diffuse interface region
* ``Ai`` -- surface area for boundary of diffuse interface region
* ``Hi`` -- integral mean curvature for boundary of diffuse interface region
* ``Xi`` -- Euler characteristic for diffuse interface region
* ``Vic`` -- volume for connected interface region
* ``Aic`` -- surface area for boundary of connected interface region
* ``Hic`` -- Integral mean curvature for connected interface region
* ``Xic`` -- Euler characteristic for connected interface region
* ``Nic`` -- number of connected components in connected interface region
