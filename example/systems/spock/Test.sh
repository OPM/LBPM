#!/bin/bash
#SBATCH -A CSC380
#SBATCH -J sphere_test
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:05:00
#SBATCH -p caar
#SBATCH -N 1

module load craype-accel-amd-gfx908
module load PrgEnv-cray
#module load rocm
module load rocm/4.2.0

export LBPM_DIR=/ccs/proj/csc380/mcclurej/spock/install/lbpm/tests
#export MPICH_RDMA_ENABLED_CUDA=1
#export MPICH_ENV_DISPLAY=1
#export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_NO_ASYNC_MEMCPY=0
export MPICH_SMP_SINGLE_COPY_MODE=CMA
#export MPICH_DBG_FILENAME="./mpich-dbg.log"
export MPICH_DBG_CLASS=ALL    
export MPICH_DBG_LEVEL=VERBOSE
export MPICH_DBG=yes
#export PMI_DEBUG=1
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
#srun -n1 --ntasks-per-node=1 --accel-bind=g --gpus-per-task=1 $LBPM_DIR/lbpm_color_simulator spheres322.db


srun -n1 --ntasks-per-node=1 --accel-bind=g --gpus-per-task=1 --verbose --export=ALL $LBPM_DIR/TestCommD3Q19 test.db

