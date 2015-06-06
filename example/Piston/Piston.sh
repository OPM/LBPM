#!/bin/bash

# This runscript is for the Virginia Tech supercomputer HokieSpeed
#PBS -l walltime=72:00:00
# Set the number of nodes, and the number of processors per node (generally should be 6)
#PBS -l nodes=1:ppn=12
##PBS -l nodes=hs195
#PBS -A arcadm
# Access group, queue, and accounting project
#PBS -W group_list=arcadm
#PBS -q large_q

module purge
module load gcc cuda mvapich2/1.9rc1

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE > hostfile

echo "------------------------------------------"
echo "Running LBM using MPI!" 
echo "Number of processors = " $PBS_NP
echo "------------------------------------------"

mpirun -np 4 ~/install-LBPM-WIA/bin/lb2_Color_wia_mpi >

exit;
