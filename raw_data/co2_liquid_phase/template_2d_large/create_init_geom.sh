#!/bin/bash

nbath=1296
nph=324

nat=$(echo $(($nbath*108 + $nph)))
echo "$nbath CO2 bath + $nph photons, in total $nat atoms"

filename="init.xyz"

rm $filename

echo "$nat" >> $filename
echo " " >> $filename
for (( c=1; c<=$nbath; c++ ))
do
echo "co2 $c"
cat co2_1b.xyz | tail -n108 >> $filename
done

for (( c=1; c<=$nph; c++ ))
do
echo "ph $c"
echo "L 0.0 0.0 0.0" >> $filename
done
