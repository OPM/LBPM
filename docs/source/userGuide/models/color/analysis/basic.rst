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
* ``krn`` -- non-wetting fluid effective permeability
* ``krn`` -- non-wetting fluid effective permeability
* ``krwf`` -- non-wetting fluid effective permeability
* ``krnf`` -- non-wetting fluid effective permeability

  
More comprehensive analysis is performed in the ``subphase`` analysis module. 

.. code-block:: bash

Analysis {
    analysis_interval = 1000
    visualization_interval = 20000
 }

