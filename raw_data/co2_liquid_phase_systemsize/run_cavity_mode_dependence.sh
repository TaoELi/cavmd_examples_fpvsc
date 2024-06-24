#!/bin/bash

# we will perform additional simulations to check the cavity mode dependence by increasing the number of cavity photon modes, leaving the size of the molecular subsystems unchanged

Nbath=36

for M in 12 16 20 24 30 #1 2 3 4 5 6 7 8 9 10
do

echo "increasing the cavity mode number by a factor of $M, system size by a factor of 4"
N=$(($Nbath * 4))
Nmode=$(($Nbath * $M))
echo "the number of molecular baths is $N"
echo "the number of cavity modes is $Nmode"

dir=mode_times_$M

mkdir $dir
cd $dir

# 1. prepare the init.xyz geometry
cp ../template/co2.xyz .
cp ../template/L.xyz .
Nat=$(($Nmode + 4*108*36))
echo "There are $Nat atoms + photons in the system"
echo "preparing the init.xyz file for the simulation"

rm init.xyz
touch init.xyz
printf "$Nat\n\n" >> init.xyz

for i in $(seq 1 $N)
do
cat co2.xyz >> init.xyz
done

for i in $(seq 1 $M)
do
cat L.xyz >> init.xyz
done

# 2. we now find tune input xml files
E0=$(echo "print(5e-5/4**0.5)" | python3)
echo "the light-matter coupling per molecule is $E0"

domega=$(echo "print(50.0/$M)" | python3)
echo "the spacing between photon in-plane frequency is $domega"

echo "the number of photonic modes is $Nmode"

# the tricky one: generate 1D array molecular grid points
x_grid_1d=" "
for i in $(seq 1 $N)
do
x=$(echo "print($i/(1.0 + $N))" | python3)

if  [ $i -eq $N ]; then
x_grid_1d="$x_grid_1d $x"
else
x_grid_1d="$x_grid_1d $x,"
fi

done

echo "$x_grid_1d"

# now replacing the parameters in i-pi inputs
cp ../template/input_eq.xml .
cp ../template/input_traj.xml.bak .

sed -i "s/0e-4/$E0/" input_eq.xml
sed -i "s/50.0 </$domega </" input_eq.xml
sed -i "s/bath> 36 </bath> 144 </" input_eq.xml
sed -i "s/x> 36 </x> $Nmode </" input_eq.xml
sed -i "s/<x_grid_1d>/<x_grid_1d> [$x_grid_1d]/" input_eq.xml

sed -i "s/0e-4/$E0/" input_traj.xml.bak
sed -i "s/50.0 </$domega </" input_traj.xml.bak
sed -i "s/bath> 36 </bath> 144 </" input_traj.xml.bak
sed -i "s/x> 36 </x> $Nmode </" input_traj.xml.bak
sed -i "s/<x_grid_1d>/<x_grid_1d> [$x_grid_1d]/" input_traj.xml.bak

# copy the other files here
cp ../template/data.lmp ../template/in.lmp ../template/submit.sh  .
sed -i "s/E0=0e-4/E0=$E0/" submit.sh
sed -i "s/RUN=size/RUN=phsize$M/" submit.sh
sed -i "s/fpcavmd/fpmode_$M/" submit.sh


# finally, submit the jobs
sbatch submit.sh

cd ..

done

