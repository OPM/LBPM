#!/bin/bash
#PBS -A GEO106
#PBS -N Drain
#PBS -j oe
#PBS -l walltime=12:00:00,nodes=432
##PBS -l gres=widow2%widow3
##PBS -q killable
##PBS -q debug

#cd /tmp/work/$USER
date

cd $PBS_O_WORKDIR


#LBPM_WIA_INSTALL_DIR=/lustre/atlas/scratch/mcclurej/geo106/FAST/LBPM-WIA
LBPM_WIA_INSTALL_DIR=/lustre/atlas/scratch/mcclurej/geo106/install-LBPM-WIA

source $MODULESHOME/init/bash
module unload PrgEnv-gnu PrgEnv-pgi PrgEnv-cray PrgEnv-intel
module load PrgEnv-gnu
module load cudatoolkit

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_MAX_THREAD_SAFETY=multiple

#LBPM_WIA_INSTALL_DIR=/lustre/atlas/scratch/mcclurej/geo106/install-LBPM-WIA

#echo "PBS_O_WORKDIR: `echo $PBS_O_WORKDIR`"
#source $MODULESHOME/init/bash
#module swap cray-mpich2 cray-mpich2/5.6.3
#module load cudatoolkit

#export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
#export MPICH_RDMA_ENABLED_CUDA=1

#LIST=$(ls | grep "pt")
LIST="pt10 pt18"

for dir in $LIST; do
    cd $dir
#    rm *.tcat
#    aprun -n 216 -N 1 $LBPM_WIA_INSTALL_DIR/bin/lbpm_sphere_pp &
    aprun -n 216 -N 1 -d 16 $LBPM_WIA_INSTALL_DIR/bin/lbpm_color_simulator &
    cd ..
    pid="$pid "$!
done

echo "Waiting for all jobs to complete"
for i in $pid; do
    wait $i
done

exit;


#    aprun -n 512 -N 1 -d 16 $LBPM_WIA_INSTALL_DIR/bin/lbpm_color_simulator >> drainage.log 
