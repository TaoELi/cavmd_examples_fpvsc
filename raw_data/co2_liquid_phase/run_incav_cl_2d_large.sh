#!/bin/bash

freqcav=2320

for E0 in 8.3e-6 #0e-4 2e-4 4e-4 6e-4
do

dir=E0_$E0\_2d

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template_2d_large/data.lmp ../template_2d_large/in.lmp ../template_2d_large/init.xyz ../template_2d_large/input_eq.xml ../template_2d_large/input_traj.xml.bak ../template_2d_large/submit.sh  .

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<E0> 0e-4/<E0> $E0/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=eq/RUN=eq-2d/" submit.sh

sbatch submit.sh

cd ..

done
