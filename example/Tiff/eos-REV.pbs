#!/bin/bash
#PBS -A GEO106
#PBS -N Juanes
#PBS -j oe
#PBS -l walltime=02:00:00,nodes=1
##PBS -l walltime=01:00:00,nodes=18
##PBS -l gres=widow2%widow3
##PBS -q killable
##PBS -q debug

#cd /tmp/work/$USER
date

cd $PBS_O_WORKDIR

LBPM_WIA_INSTALL_DIR=/ccs/proj/geo106/eos/LBPM-WIA

source $MODULESHOME/init/bash
module swap PrgEnv-intel PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load python_anaconda

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

## This is the code used to extract the tifs to subdirectories
for i in `seq -w 0 804`; do 
    file="ConcatenatedStacks0"$i".tif"; cp $file SubSet1; 
done

for i in `seq -w 805 1609`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet2; 
done

for i in `seq -w 1610 2414`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet3; 
done

for i in `seq -w 2415 3219`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet4; 
done

for i in `seq -w 3220 4024`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet5; 
done

for i in `seq -w 4025 4829`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet6; 
done

for i in `seq -w 4830 5634`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet7; 
done

for i in `seq -w 5635 6439`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet8; 
done

for i in `seq -w 6440 7244`; do 
    file="ConcatenatedStacks"$i".tif"; cp $file SubSet9; 
done


# This is the code used to convert tif to raw 8-bit binary
LIST=$(ls | grep SubSet)

for dir in $LIST; do 
   echo "Processing $dir"
   cd $dir
   python $LBPM_WIA_INSTALL_DIR/example/Tiff/TiffPreprocessor.py
   cd ..
done


#     aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_decomp 0 2

#     aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_segmented_pp

exit;
