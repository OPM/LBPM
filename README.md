LBPM-WIA
========

Lattice Boltzmann Methods for Porous Media with Integrated Averaging

========

Notes on installation

* example configure scripts for cmake are in the sample_scripts directory
* required dependencies - MPI, C++ 11
* optional dependencies - NetCDF, CUDA, TimerUtility

Configure, build & install procedure
1) create an empty directory to install (do not build in source!)
   mkdir /path/to/my/install
   cd /path/to/my/install
2) edit configure script from sample_scripts directory and configure (e.g.)
   /path/to/LBPM-WIA/sample_scripts/configure_desktop 
3) compile and install
   make -j4 && make install
4) run the unit tests to make sure it works
   ctest

