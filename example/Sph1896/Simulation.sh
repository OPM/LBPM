#!/bin/bash

LBPM_INSTALL_DIR=../../bin

mpirun -np 8 $LBPM_INSTALL_DIR/GenerateSphereTest input.db
mpirun -np 8 $LBPM_INSTALL_DIR/lbpm_color_simulator input.db
