#!/bin/bash

freqcav=2320

for E0 in 5e-5 #0e-4 2e-4 4e-4 6e-4
do

dir=E0_$E0\_excited2

echo "working at $dir"

mkdir $dir

cd $dir

cp ../template_excited/data.lmp ../template_excited/in.lmp ../template_excited/input_traj.xml.bak ../template_excited/submit.sh  .

for i in {1..40}
do
   tail -n 3926 ../E0_$E0/simu_$i.xc.xyz > init_$i\.xyz
   sed -i '3902s/.*/       L  0.0000000e-01 1.0000000e+02  0.000000e-01/' init_$i\.xyz
done

# modify light-matter coupling
sed -i "s/<E0> 0e-4/<E0> $E0/" input_traj.xml.bak
sed -i "s/<omega_c_cminv> 2320.0/<omega_c_cminv> $freqcav/" input_traj.xml.bak
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=eq/RUN=excited2/" submit.sh


sbatch submit.sh

cd ..

done
