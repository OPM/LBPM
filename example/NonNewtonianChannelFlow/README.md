# 1. Edit the Domain.in file based on the desired system size


# 2. Run the pre-processor 
export TUBEWIDTH=20
export LAYERWIDTH=0
mpirun -np 1 ../../tests/lbpm_plates_pp $TUBEWIDTH $LAYERWIDTH


# 3. Run single-phase simulation within the domain
mpirun -np 1 ../../tests/lbpm_permeability_simulator
