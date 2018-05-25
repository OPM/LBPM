#!/bin/bash

# Run morphological analysis and set up input files for media
LABEL="sph1964"
#aprun -n 64 $LBPM_WIA_INSTALL_DIR/bin/lbpm_morph_pp $r > morph.log
echo "Quartile Radius" > poresize.quartiles
grep -A4 Quartile morph.log | tail -4 >> poresize.quartiles
python ~/LBPM-WIA/workflows/imbibition/setup-imbibition.py $LABEL

FILES=$(ls | grep "drain_")

for src in $FILES; do  
     dest=$(echo $src | sed 's/drain_/imb_/g')
     echo "Initializing $dest from primary drainage"
     cp -r $src $dest
     cp Color.in.imbibition $dest/Color.in
done

exit;

