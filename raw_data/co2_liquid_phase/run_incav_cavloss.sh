#!/bin/bash

freqcav=2320

E0=5e-5

for cavlifetime in 300
do

dir=cavlifetime_$cavlifetime

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template/data.lmp ../template/in.lmp ../template/init.xyz ../template/input_eq.xml ../template/input_traj_cavloss.xml.bak ../template/submit.sh  .

mv input_traj_cavloss.xml.bak input_traj.xml.bak

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<E0> 0e-4/<E0> $E0/" input_eq.xml
sed -i "s/<omega_c units='inversecm'> 2320.0/<omega_c units='inversecm'> $freqcav/" input_eq.xml
sed -i "s/<omega_c units='inversecm'> 2320.0/<omega_c units='inversecm'> $freqcav/" input_traj.xml.bak
sed -i "s/<tau_l units='femtosecond'> 100.0/<tau_l units='femtosecond'> $cavlifetime/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=eq/RUN=loss$cavlifetime/" submit.sh

sbatch submit.sh

cd ..

done
