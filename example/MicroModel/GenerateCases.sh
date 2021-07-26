#!/bin/bash


tau1=1.18
tau2=0.7
alpha=0.005
Q="1.179757e-08 1.179757e-07 1.179757e-06"

#Cases for drainage
DrainWet="0.79 0.47 0.0" 

# Cases for imbibition
ImbWet="0.92 0.47" 

for q in $Q; do 
  echo $q;
  flux=$(echo $q | sed 's/1.179757e-08/0.002/g')
  flux=$(echo $flux | sed 's/1.179757e-07/0.02/g')
  flux=$(echo $flux | sed 's/1.179757e-06/0.2/g')

   for i in $DrainWet; do
       NAME="Juanes_drain_Q"$flux"_wet"$i
       echo $NAME
       mkdir $NAME
       echo "$tau1 $tau2" > $NAME/Color.in
       echo "$alpha 0.95 $i" >> $NAME/Color.in
       echo "0.0" >> $NAME/Color.in
       echo "0.0 0.0 0.0" >> $NAME/Color.in
       echo "0 4 $q 0.0" >> $NAME/Color.in
       echo "5000000 25000 1e-5" >> $NAME/Color.in
       echo "1000" >> $NAME/Color.in
   done

   for i in $ImbWet; do
       NAME="Juanes_imb_Q"$flux"_wet"$i
       echo $NAME
       mkdir $NAME
       echo "$tau1 $tau2" > $NAME/Color.in
       echo "$alpha 0.95 $i" >> $NAME/Color.in
       echo "0.0" >> $NAME/Color.in
       echo "0.0 0.0 0.0" >> $NAME/Color.in
       echo "0 4 $q 0.0" >> $NAME/Color.in
       echo "5000000 25000 1e-5" >> $NAME/Color.in
       echo "1000" >> $NAME/Color.in
   done
  
done