#!/bin/bash

#SBATCH -A CSC380
#SBATCH -J MPI-singlenode
#SBATCH -o %x-%j.out
#SBATCH -t 0:10:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --exclusive

# MODULE ENVIRONMENT
module load PrgEnv-amd
module load rocm/4.5.0
module load cray-mpich
module load cray-hdf5-parallel
#module load craype-accel-amd-gfx908

## These must be set before compiling so the executable picks up GTL
export PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
export PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
export MPICH_GPU_SUPPORT_ENABLED=1

#export MPL_MBX_SIZE=1024000000

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

export LBPM_BIN=/ccs/proj/csc380/mcclurej/crusher/LBPM/tests

echo "Running Color LBM"

MYCPUBIND="--cpu-bind=verbose,map_cpu:57,33,25,1,9,17,41,49"

srun --verbose -N1 -n8 --cpus-per-gpu=8 --gpus-per-task=1 --gpu-bind=closest ${MYCPUBIND} $LBPM_BIN/TestCommD3Q19 multinode.db
#srun --verbose -N1 -n2 --mem-per-gpu=8g --cpus-per-gpu=1 --gpus-per-node=2 --gpu-bind=closest $LBPM_BIN/lbpm_permeability_simulator input.db

exit;
