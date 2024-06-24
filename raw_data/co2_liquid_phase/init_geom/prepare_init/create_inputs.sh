#!/bin/bash

while getopts n:f:w: option
do
case "${option}"
in
n) NMOLECULE=${OPTARG};;
f) LAMMPS_DATA_FILE=${OPTARG};;
w) SIMU_CELL_LENGTH=${OPTARG};;
esac
done

NMOLECULE=36
LAMMPS_DATA_FILE=data_36co2.lmp
SIMU_CELL_LENGTH=25.263

cd generate_initial_config

./generate_initial_config.sh -n $NMOLECULE -f $LAMMPS_DATA_FILE -w $SIMU_CELL_LENGTH

cd ..

cd generate_lammps_inputs

./generate_lmp_data.sh -n $NMOLECULE -f $LAMMPS_DATA_FILE -w $SIMU_CELL_LENGTH

cd ..

cd minimize_initial_config

./generate_opt_geometry.sh -n $NMOLECULE -f $LAMMPS_DATA_FILE -w $SIMU_CELL_LENGTH

cd ..

# We do a final modification
mv $LAMMPS_DATA_FILE data.lmp
mv in_opt.lmp in.lmp
sed -i "s/$LAMMPS_DATA_FILE/data.lmp/" in.lmp

#sed -i "s/343/$NMOLECULE/" create_photon_params.py
#python create_photon_params.py
#cp test.json photon_params.json

#sed -i "s/53.557/$SIMU_CELL_LENGTH/g" input_equlibrate.xml
