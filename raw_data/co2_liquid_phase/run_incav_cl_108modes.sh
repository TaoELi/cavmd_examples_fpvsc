#!/bin/bash

freqcav=2320

for E0 in 5e-5 #0e-4 2e-4 4e-4 6e-4
do

dir=E0_$E0\_108modes

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template/data.lmp ../template/in.lmp ../template/init.xyz ../template/input_eq.xml ../template/input_traj.xml.bak ../template/submit.sh  .

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<E0> 0e-4/<E0> $E0/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_traj.xml.bak
sed -i "s/<domega_x_cminv> 50.0/<domega_x_cminv> 16.67/" input_eq.xml
sed -i "s/<domega_x_cminv> 50.0/<domega_x_cminv> 16.67/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=eq/RUN=eq108/" submit.sh
sed -i "s/pi_hammes_schiffer/week/" submit.sh
sed -i "1s/3924/3996/" init.xyz
sed -i "s/<n_mode_x> 36/<n_mode_x> 108/" input_eq.xml
sed -i "s/<n_mode_x> 36/<n_mode_x> 108/" input_traj.xml.bak

sbatch submit.sh

cd ..

done
