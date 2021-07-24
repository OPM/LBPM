======================================
Color model -- Basic Analysis
======================================

The basic analysis routine for the LBPM color model logs a time series
of averaged information to the space-delimited CSV file ``timelog.csv``.
The ``analysis_interval`` key should be specified to control the interval at
to perform basic analysis, which corresponds to the number of simulation timesteps
to perform between analyses. The basic analysis routine is designed to
be lightweight so that it can be performed frequently, and it is almost always
useful to enable it. Results  will be immediately logged
to ``timelog.csv`` which can be used to assess whether or not the simulation is
behaving as expected. Furthermore, the analysis framework will check
simulation measures to verfiy that no ``NaN`` are detected for the fluid
pressure and flow rate.

The list of measures logged to ``timelog.csv`` are defined as follows.
The region of space occupied by the wetting fluid is determined from the
phase indicator field, :math:`\Omega_w:\phi<0` 


* ``sw`` -- water saturation (fluid component 2)
* ``krw`` -- water effective permeability 
* ``krw`` -- effective permeability for fluid w
* ``krn`` -- effective permeability for fluid n
* ``krwf`` -- effective permeability for fluid w (with film uncertainty estimate)
* ``krnf`` -- effective permeability for fluid n (with film uncertainty estimate)
* ``vw`` -- speed for the non-wetting fluid
* ``vn`` -- speed for the wetting fluid
* ``force`` -- magnitude for effective driving force
* ``pw`` -- average pressure for fluid w
* ``pn`` -- average pressure for fluid n
* ``wet`` -- total solid wetting energy

  
More comprehensive analysis is performed in the ``subphase`` analysis module. 


