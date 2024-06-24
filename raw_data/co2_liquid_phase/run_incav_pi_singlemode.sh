#!/bin/bash

freqcav=2320

for E0 in 5e-5 #0e-4 2e-4 4e-4 6e-4
do

dir=pi_E0_$E0\_singlemode

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template_pi/data.lmp ../template_pi/in.lmp ../template_pi/init.xyz ../template_pi/input_eq.xml ../template_pi/input_traj.xml.bak ../template_pi/submit.sh  .

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<E0> 0e-4/<E0> $E0/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=eq/RUN=eq-pisingle/" submit.sh
sed -i "1s/3924/3889/" init.xyz
sed -i "s/<n_mode_x> 36/<n_mode_x> 1/" input_eq.xml
sed -i "s/<n_mode_x> 36/<n_mode_x> 1/" input_traj.xml.bak

sbatch submit.sh

cd ..

done
