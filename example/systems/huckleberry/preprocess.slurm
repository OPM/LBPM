#!/bin/bash
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 0:30:00
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH --account=slurmtest

echo "LBPM on huckleberry"

module load gcc/7.3.0 openmpi/3.1.2 cuda 
module load szip 
module load zlib 
module list
module load hdf5/1.8.12 
module load silo 
#module load lbpm

export LBPM_BIN=$HOME/install/huckleberry/LBPM/bin
export CUDA_VISIBLE_DEVICES=0,1,2,3

export OMPI_MCA_btl_smcuda_use_cuda_ipc=0
export OMPI_MCA_pml=ob1
export OMPI_MCA_mpi_warn_on_fork=0

#srun --exclusive -n 4 --cpu-bind=v,map_cpu=0,40,80,120 --accel-bind=g --gres=gpu:4 $LBPM_BIN/lbpm_color_simulator input.db
#
#

mpirun -np 1 $LBPM_BIN/lbpm_serial_decomp input.db
mpirun --oversubscribe --npernode 8 -mca pml ob1 -mca btl_smcuda_use_cuda_ipc 0 --bind-to core:overload-allowed  $LBPM_BIN/lbpm_morphdrain_pp input.db

exit;
