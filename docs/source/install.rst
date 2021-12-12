============
Installation
============

Sample scripts to configure and build LBPM are included in the `sample_scripts` directory. To build this package, MPI and cmake are required, along with a compiler that supports the C++-14 standard. Building in source is not supported.

The essential dependencies needed to build LBPM are:

1.  cmake (version 3.9 or higher)
2.  MPI
3.  HDF5
4.  silo

*************************
Building Dependencies 
*************************
(skip if they are already installed on your system)

1. Download third-party library source files

.. code-block:: bash

   zlib-1.2.11.tar.gz
   hdf5-1.8.12.tar.gz 
   silo-4.10.2.tar.gz 


2. Set up path where you want to install 

.. code-block:: bash

   export MPI_DIR=/path/to/mpi
   export LBPM_ZLIB_DIR=/path/to/zlib
   export LBPM_HDF5_DIR=/path/to/hdf5
   export LBPM_SILO_DIR=/path/to/silo

3. Build third-party library dependencies 

.. code-block:: bash

   tar -xzvf zlib-1.2.11.tar.gz 
   tar -xzvf hdf5-1.8.12.tar.gz 
   tar -xzvf silo-4.10.2.tar.gz 
   cd zlib-1.2.11

   ./configure --prefix=$LBPM_ZLIB_DIR && make && make install

   cd ../hdf5-1.8.12

   CC=$MPI_DIR/bin/mpicc  CXX=$MPI_DIR/bin/mpicxx CXXFLAGS="-fPIC -O3 -std=c++14" \
   ./configure --prefix=$LBPM_HDF5_DIR --enable-parallel --enable-shared --with-zlib=$LBPM_ZLIB_DIR 
    make && make install


   cd ../silo-4.10.2

   CC=$MPI_DIR/bin/mpicc  CXX=$MPI_DIR/bin/mpicxx CXXFLAGS="-fPIC -O3 -std=c++14" \
   ./configure --prefix=$LBPM_SILO_DIR -with-hdf5=$LBPM_HDF5_DIR/include,$LBPM_HDF5_DIR/lib --enable-static 
    make && make install
 
*************************
Building LBPM 
*************************

Many HPC systems will include all of these dependencies, and LBPM can be built simply by setting the paths to the associated libraries in the cmake configure line. 

The steps typically used to build LBPM are as follows:

1. Set environment variables for the source directory `$LBPM_WIA_SOURCE` and build directory `$LBPM_WIA_DIR`  

.. code-block:: bash

   export LBPM_SOURCE=/path/to/source/LBPM  
   export LBPM_DIR=/path/to/build/LBPM

2. Set environment variables for the path to HDF5 and SILO (required),  and optionally for TimerUtility and NetCDF (optional)  

.. code-block:: bash

   export LBPM_HDF5_DIR=/path/to/hdf5  
   export LBPM_SILO_DIR=/path/to/silo 
   export LBPM_TIMER_DIR=/path/to/timer 
   export LBPM_NETCDF_DIR=/path/to/netcdf

3. Create the build directory and navigate to it  

.. code-block:: bash

   mkdir $LBPM_WIA_DIR  
   cd $LBPM_WIA_DIR

4. Configure the project. Numerous scripts exist to build LBPM on different HPC clusters, which are available in the  `$LBPM_SOURCE/sample_scripts/` directory. It is often possible to run these scripts directly if one exists for a system similar to the one you are building on. For a standard CPU build:  

.. code-block:: bash

    cmake                                           \  
        -D CMAKE_BUILD_TYPE:STRING=Release          \   
        -D CMAKE_C_COMPILER:PATH=mpicc              \   
        -D CMAKE_CXX_COMPILER:PATH=mpicxx           \  
        -D CMAKE_C_FLAGS="-fPIC"                    \  
        -D CMAKE_CXX_FLAGS="-fPIC"                  \  
        -D CMAKE_CXX_STD=14                         \  
        -D USE_TIMER=0                              \  
            -D TIMER_DIRECTORY=$LBPM_TIMER_DIR     \  
        -D USE_NETCDF=0                             \  
            -D NETCDF_DIRECTORY=$LBPM_NETCDF_DIR   \  
        -D USE_SILO=1                               \  
           -D HDF5_DIRECTORY=$LBPM_HDF5_DIR         \  
           -D SILO_DIRECTORY=$LBPM_SILO_DIR         \  
        -D USE_CUDA=0                               \
        $LBPM_SOURCE  

For GPU support, it is necessary to have CUDA along with a GPU-aware MPI implementation. Otherwise, the LBPM routines should behave identically irrespective of the underlying hardware.  

.. code-block:: bash

    cmake                                           \  
        -D CMAKE_BUILD_TYPE:STRING=Release          \   
        -D CMAKE_C_COMPILER:PATH=mpicc              \   
        -D CMAKE_CXX_COMPILER:PATH=mpicxx           \  
        -D CMAKE_C_FLAGS="-fPIC"                    \  
        -D CMAKE_CXX_FLAGS="-fPIC"                  \  
        -D CMAKE_CXX_STD=14                         \  
        -D USE_TIMER=0                              \  
            -D TIMER_DIRECTORY=$LBPM_TIMER_DIR     \  
        -D USE_NETCDF=0                             \  
            -D NETCDF_DIRECTORY=$LBPM_NETCDF_DIR   \  
        -D USE_SILO=1                               \  
           -D HDF5_DIRECTORY=$LBPM_HDF5_DIR         \  
           -D SILO_DIRECTORY=$LBPM_SILO_DIR         \  
        -D USE_CUDA=1                               \  
        -D CMAKE_CUDA_FLAGS="-arch sm_70"           \  
        $LBPM_SOURCE   

5. Build the project (using four cores to build)   

.. code-block:: bash
  
   make -j4

6. Install the project  

.. code-block:: bash

   make install

7. Run the tests to make sure they execute correctly (on a cluster, it is recommended to run these using the batch system rather than on the head node)  

.. code-block:: bash

   ctest


*************************
Sample Scripts
*************************

The LBPM repository contains sample scripts showing successful CMake configuration, build and
install steps for a range of systems. Refer to the project sub-directory below for these examples.

.. code-block:: bash

    ls $LBPM_SOURCE/sample_scripts
    
