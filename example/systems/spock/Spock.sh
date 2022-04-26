#!/bin/bash
#SBATCH -A CSC380
#SBATCH -J sphere_test
#SBATCH -o %x-%j.out
#SBATCH -t 00:05:00
#SBATCH -p caar
#SBATCH -N 1

module load rocm/4.2.0
export LBPM_DIR=/ccs/proj/csc380/mcclurej/spock/install/lbpm/tests
export MPICH_SMP_SINGLE_COPY_MODE=CMA

#srun -n1 --ntasks-per-node=1 --accel-bind=g --gpus-per-task=1 $LBPM_DIR/lbpm_color_simulator spheres322.db


srun -n1 --ntasks-per-node=1 --accel-bind=g --gpus-per-task=1 $LBPM_DIR/TestCommD3Q19 spheres322.db

