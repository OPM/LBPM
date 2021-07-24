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

Since it is more computationally expensive to carry out these operations compared to the
basic analysis, it may be useful to choose ``subphase_analysis_interval`` to be larger than
``analysis_interval``. Nevertheless, since analysis is performed in memory, it is orders of
magnitude faster than to analyze data in situ rather than writing fields to disc. Monitoring
the performance in MLUPS provides an easy way to tune the analysis interval so that the
overall simulation performance is not subject to significant penalties. 

There are three main entities that are constructed for subphase analysis

* :math:`\Omega_w:\phi<0` : region of space filled by fluid w
* :math:`\Omega_n:\phi<0` : region of space filled by fluid n
* :math:`\Omega_i: |\nabla \phi<0| > \epsilon` : region of space for the interface region

The phase regions are defined identical to what is reported in ``timelog.csv``.
The interface region is defined explicitly as the portion of space where
significant composition gradients are present. 

The list of measures logged to ``subphase.csv`` are defined as follows.
The region of space occupied by the wetting fluid is determined from the
phase indicator field, :math:`\Omega_w:\phi<0` 

* ``sw`` -- water saturation (fluid component 2)
* ``krw`` -- water effective permeability 
* ``krn`` -- non-wetting fluid effective permeability
* ``krn`` -- non-wetting fluid effective permeability
* ``krwf`` -- non-wetting fluid effective permeability
* ``krnf`` -- non-wetting fluid effective permeability

