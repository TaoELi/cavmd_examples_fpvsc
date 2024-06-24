#!/bin/bash

# we will perform additional simulations to check the system size dependence by increasing the number of cavity photon modes and molecular subsystems at the same time

Nbath=36

for M in 1 2 4 6 8
do

echo "increasing the cavity mode size by a factor of $M"

dir_eq=mode_times_$M
dir=excite_ph_mode_times_$M

mkdir $dir
cd $dir

# 1. prepare the initital geometry by displacing the y-coordinate of the No.12 photonic coordinate
for checkpoint in ../$dir_eq/init*checkpoint
do
    echo "changing initial condition for $checkpoint"
    python ../change_initial_condition_excite_ph.py $checkpoint 4 $M
done

# 2. copy the other files to here
cp ../$dir_eq/input_traj.xml.bak .
cp ../$dir_eq/in.lmp ../$dir_eq/data.lmp .
cp ../template/submit_noneq.sh .

# add the output of input_traj.xml.bak
awk 'NR==5{print "<trajectory filename=\"vc\" stride=\"4\" format=\"xyz\"> v_centroid_ph </trajectory>"}1' input_traj.xml.bak > hehe
mv hehe input_traj.xml.bak

E0=$(echo "print(5e-5/$M**0.5)" | python3)
echo "the light-matter coupling per molecule is $E0"
sed -i "s/E0=0e-4/E0=$E0/" submit_noneq.sh
sed -i "s/RUN=ex/RUN=exm/" submit_noneq.sh

# finally, submit the jobs
sbatch submit_noneq.sh

cd ..

done
