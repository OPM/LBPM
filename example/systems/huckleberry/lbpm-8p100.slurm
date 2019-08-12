#!/bin/bash
#SBATCH -p normal_q
#SBATCH -N 2
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH --mem=500G
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH --account=slurmtest

echo "LBPM on huckleberry"

module load gcc/7.3.0 openmpi/3.1.2 cuda 
module load szip 
module load zlib 
module load hdf5/1.8.12 
module load silo 
#module load lbpm
#module list

export LBPM_BIN=$HOME/install/huckleberry/LBPM/bin

export CUDA_VISIBLE_DEVICES=0,1,2,3
export OMPI_MCA_btl_smcuda_use_cuda_ipc=0
export OMPI_MCA_btl_mpi_common_cuda_verbose=100
export OMPI_MCA_pml=ob1
export OMPI_MCA_mpi_warn_on_fork=0

rm host.list
for gpu in `seq 1 4`; do
  scontrol show hostname $SLURM_NODELIST >> host.list
done

MPIARGS="-x OMPI_MCA_pml --mca mpi_common_cuda_verbose 100 --mca pml ob1 --mca btl_smcuda_use_cuda_ipc 0 --bind-to core --hostfile host.list"


MPIARGS="--hostfile host.list --bind-to core --mca pml ob1 --mca btl vader,self,smcuda,openib  --mca btl_openib_warn_default_gid_prefix 0  --mca btl_smcuda_use_cuda_ipc_same_gpu 0  --mca btl_openib_want_cuda_gdr 0  --mca btl_openib_cuda_async_recv false --mca btl_smcuda_use_cuda_ipc 0 --mca btl_openib_allow_ib true --mca btl_openib_cuda_rdma_limit 1000 -x LD_LIBRARY_PATH"
#mpirun -np 8 -mca pml ob1 -mca btl_smcuda_use_cuda_ipc 0 --bind-to core:overload-allowed $LBPM_BIN/lbpm_color_simulator input.db
#mpirun -np 8 $MPIARGS  $LBPM_BIN/lbpm_morphdrain_pp input.db
mpirun -np 8 $MPIARGS  $LBPM_BIN/lbpm_color_simulator input.db

exit;
