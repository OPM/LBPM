#!/bin/bash
#PBS -A GEO106
#PBS -N Random
#PBS -j oe
##PBS -l walltime=02:00:00,nodes=216
#PBS -l walltime=12:00:00,nodes=18
##PBS -l gres=widow2%widow3
##PBS -q killable
##PBS -q debug

#cd /tmp/work/$USER
date

cd $PBS_O_WORKDIR

#LBPM_WIA_INSTALL_DIR=/lustre/atlas/proj-shared/geo106/build-eos-LBPM-WIA

LBPM_WIA_INSTALL_DIR=/ccs/proj/geo106/eos/LBPM-WIA

#echo "PBS_O_WORKDIR: `echo $PBS_O_WORKDIR`"
source $MODULESHOME/init/bash
module swap PrgEnv-intel PrgEnv-gnu
module load python_anaconda

#module swap cray-mpich2 cray-mpich2/5.6.3
#module load cudatoolkit

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

LABEL=$(basename $PWD)

NUMPROCS=288
# Generate a bunch of random states
for sat in `seq -w 8 4 92`; do 
      cp ../MEDIA/ID* ./
      sw="0.$sat"
      aprun -n $NUMPROCS $LBPM_WIA_INSTALL_DIR/bin/lbpm_random_pp $sw 1
  
      tag="sw_$sw"
      DIR=$LABEL"_random_"$tag
      mkdir -p $DIR
      cp ID.* $DIR
      cp SignDist.* $DIR
      cp Domain.in $DIR
      cp Color.in $DIR
done

exit;

