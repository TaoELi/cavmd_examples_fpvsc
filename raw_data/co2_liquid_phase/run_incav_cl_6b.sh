#!/bin/bash

freqcav=2320

for E0 in 1.225e-4 #0e-4 2e-4 4e-4 6e-4
do

dir=E0_$E0\_6b

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template_6baths/data.lmp ../template_6baths/in.lmp ../template_6baths/init.xyz ../template_6baths/input_eq.xml ../template_6baths/input_traj.xml.bak ../template_6baths/submit.sh  .

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<E0> 0e-4/<E0> $E0/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_eq.xml
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh

sbatch submit.sh

cd ..

done
