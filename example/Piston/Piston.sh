#!/bin/bash

# Lines assigning various pressure BC for Color.in
echo "0 1 1.01 0.99" > Color.in.pressures
echo "0 1 1.0125 0.9875" >> Color.in.pressures
echo "0 1 1.015 0.985" >> Color.in.pressures
echo "0 1 1.02 0.98" >> Color.in.pressures
echo "0 1 1.025 0.975" >> Color.in.pressures
echo "0 1 1.03 0.97" >> Color.in.pressures

for i in `seq 1 6`; do 
    # Set up cases for each boundary pressure pair
    dir="Case"$i
    echo $dir
    mkdir -p $dir
    # copy the domain file
    cp Domain.in $dir

    # set up each case -- parameters are fixed in Color.in, with multiple cases to set the boundary pressure
    sed -n '1p' Color.in > $dir/Color.in
    sed -n '2p' Color.in >> $dir/Color.in
    sed -n '3p' Color.in >> $dir/Color.in
    sed -n '4p' Color.in >> $dir/Color.in
#    sed -n '5p' Color.in >> $dir/Color.in
    # print the pressure values into the input file
    sed -n "${i}p" Color.in.pressures >> $dir/Color.in
    sed -n '6p' Color.in >> $dir/Color.in

done

# simulations should be run using the following syntax
# PRE-PROCESSOR - set the radius to 18 voxel lengths
#mpirun -np 10 ~/install-LBPM-WIA/bin/lbpm_captube_pp 18 1
# RUN THE SIMULAUTION 
#mpirun -np 10 ~/install-LBPM-WIA/bin/lbpm_color_simulator

exit;
