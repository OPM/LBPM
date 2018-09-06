#!/bin/bash

LBPM_DIR=../../tests

python ParallelPlates.py
mpirun -np 1 $LBPM_DIR/lbpm_serial_decomp input.db
mpirun -np 4 $LBPM_DIR/lbpm_color_simulator input.db
