#!/bin/bash

NRANKS=3600
echo $NRANKS

DIR=$NRANKS"p"
mkdir -p $DIR

BASEDIST="SignDist.0"
BASEID="ID.0"

for i in `seq -w 0 $NRANKS`; do idfile="$BASEID$i"; echo $idfile; cp ID.00000 $DIR/$idfile; done
for i in `seq -w 0 $NRANKS`; do distfile="$BASEDIST$i"; echo $distfile; cp SignDist.00000 $DIR/$distfile; done
