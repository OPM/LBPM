#!/bin/bash
#SBATCH -A CSC380
#SBATCH -J sphere_test
#SBATCH -o %x-%j.out
#SBATCH -t 00:05:00
#SBATCH -p caar
#SBATCH -N 1

module load rocm/4.2.0
export LBPM_DIR=/ccs/proj/csc380/mcclurej/spock/install/lbpm/tests

srun -n1 --ntasks-per-node=1 $LBPM_DIR/GenerateSphereTest input.db


