#!/bin/bash

#SBATCH --job-name=fpcavmd
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard

RUN=ex
E0=5e-05
    
PATTERN='<step>40000</step>'

# Here, we prepare the input files for each trajectory
for traj in {1..40}
do
    echo "Dealing with $traj slice"
    cp input_traj.xml.bak input_traj_$traj.xml
    cp in.lmp in_$traj.lmp
    # change the i-pi job name for each traj
    sed -i "s/mesitylene-pimd.1/co2.$RUN.E0_$E0.traj_$traj/" input_traj_$traj.xml
    sed -i "s/mesitylene-pimd.1/co2.$RUN.E0_$E0.traj_$traj/" in_$traj.lmp
    # change the input file for each traj 	
    sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_traj_$traj.xml
    sed -i "s/'simu'/'simu_$traj'/g" input_traj_$traj.xml

	CHECKPOINT=simu_$traj.checkpoint
	if grep -q $PATTERN $CHECKPOINT; then
	    echo "found checkpoint finished"
	    echo "Skip $traj-th sequential job"
	else
	    echo "not found checkpoint finished"
	    a=$(wc -c < $CHECKPOINT)
	    if [ ! -f "$CHECKPOINT" ] || [ $a -le 1 ]; then
	       echo "Performing $traj-th simulation 40 ps"
	       i-pi input_traj_$traj.xml &> log_$traj &
	       sleep 100s
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log
	       sleep 100s
	    else
	       echo "Continuing $traj-th simulation 20 ps"
	       i-pi $CHECKPOINT &> log_$traj &
	       sleep 100s
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log &
	       lmp < in_$traj.lmp &> info_lmp_$traj.log
	       sleep 100s
	    fi
	fi
done
wait 

