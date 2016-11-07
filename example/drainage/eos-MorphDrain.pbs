#!/bin/bash
#PBS -A GEO106
#PBS -N MorphDrain
#PBS -j oe
##PBS -l walltime=02:00:00,nodes=216
#PBS -l walltime=01:00:00,nodes=4
##PBS -l gres=widow2%widow3
##PBS -q killable
##PBS -q debug

#cd /tmp/work/$USER
date

cd $PBS_O_WORKDIR

#LBPM_WIA_INSTALL_DIR=/lustre/atlas/proj-shared/geo106/build-eos-LBPM-WIA

LBPM_WIA_INSTALL_DIR=$MEMBERWORK/geo106/eos-LBPM-WIA

#echo "PBS_O_WORKDIR: `echo $PBS_O_WORKDIR`"
source $MODULESHOME/init/bash
module swap PrgEnv-intel PrgEnv-gnu

#module swap cray-mpich2 cray-mpich2/5.6.3
#module load cudatoolkit

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
#export MPICH_RDMA_ENABLED_CUDA=1
#aprun -n 27 -N 1 /lustre/atlas/scratch/mcclurej/geo106/LBPM-WIA/bin/lb2_Color_wia_mpi

#LIST=$(ls|grep sw)
#for DIR in $LIST; do 
#  cd $DIR
#  sat=$(cat Color.in  |head -3 | tail -1)
#  aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_decomp
#  aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_pp
#  aprun -n 144 $LBPM_WIA_INSTALL_DIR/bin/lbpm_random_pp 1 0 
#  aprun -n 144 $LBPM_WIA_INSTALL_DIR/bin/lbpm_random_pp 1 1

#  aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_decomp 1 2
#  aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_pp 

#  aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_random_pp 0 1

#mkdir -p MEDIA
#cp ID* MEDIA/

cp MEDIA/ID* ./

LABEL="lrc32_sandpack"
# Run morphological drainage and set up input files for media
RADIUS="5.4 5.35 5.2 5 4.8 4.6 4.4 4.2 4 3.67 3.33 3 2.67 2.33 2 1.5"
echo "radius sw" > morphdrain.csv
for r in $RADIUS; do 
      aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_morphdrain_pp $r > morph.log
      sat=$(grep sat morph.log | sed 's/Final saturation=//g')
      tag=$(wc -l morphdrain.csv | cut -c1-2) 
      echo "$r $sat" >> morphdrain.csv
      DIR=$LABEL"_drain_"$tag
      mkdir -p $DIR
      cp ID.* $DIR
      cp SignDist.* $DIR
      cp Domain.in $DIR
done

exit;

