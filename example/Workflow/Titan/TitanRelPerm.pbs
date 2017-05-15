#!/bin/bash
#PBS -A GEO106
#PBS -N Ryan-Sandstone-RelPerm
#PBS -j oe
#PBS -l walltime=6:00:00,nodes=10660

# Number of nodes should be [# directories] x [GPU per directory]
# 

#cd /tmp/work/$USER
date

cd $PBS_O_WORKDIR

LBPM_WIA_INSTALL_DIR=/ccs/proj/geo106/titan/LBPM-WIA

source $MODULESHOME/init/bash
module unload PrgEnv-gnu PrgEnv-pgi PrgEnv-cray PrgEnv-intel
module load PrgEnv-gnu
module load cudatoolkit
module load dynamic-link
module load python wraprun

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_MAX_THREAD_SAFETY=multiple


LIST=$(ls | grep "s5_")

NUMSIM=$(ls | grep "s5_" | wc -l)

# Requested number of nodes must be 54 x $NUMSIM !!
echo "Number of simulations to run is $NUMSIM (288 GPU per simulation)" 

NUMPROCS=288

# BUIlD THE LAUNCH COMMAND 
for dir in $LIST; do
    ARGS=$ARGS" : -n $NUMPROCS -N 1 -d 16 --w-cd $dir $LBPM_WIA_INSTALL_DIR/bin/lbpm_color_simulator"
done

echo $ARGS | sed 's/./wraprun/' > launch.log

LAUNCH=$(cat launch.log)

echo "running job "
$LAUNCH

exit;


#    aprun -n 512 -N 1 -d 16 $LBPM_WIA_INSTALL_DIR/bin/lbpm_color_simulator >> drainage.log 
