#!/bin/bash

#SBATCH --job-name=fpcavmd_4
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=taoeli

RUN=size
E0=2.5e-05
    
# First, we run equilibrium simulations

FILE_EQU=init_0.checkpoint
EQCHECKPOINT=eq.checkpoint
sed -i "s/mesitylene-pimd.1/co2.$RUN.E0_$E0/" input_eq.xml
cp in.lmp in_eq.lmp
sed -i "s/mesitylene-pimd.1/co2.$RUN.E0_$E0/" in_eq.lmp
# Further modifications ...
if [ ! -f "$FILE_EQU" ]; then
    # Run for equilibrium
    # first check if the equilibrating job is not finished
    if [ -f "$EQCHECKPOINT" ]; then
        echo "Continuing equilibrating job 20 ps"
        i-pi $EQCHECKPOINT &> log_eq &
        sleep 10s
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp
        sleep 10s
    else
        echo "Performing equilibrating job 20 ps"
        i-pi input_eq.xml &> log_eq &
        sleep 10s
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp &
        lmp < in_eq.lmp
        sleep 10s
    fi
else
    echo "Skip equilibrating job"
fi

# copy the final configuration as the initial configuration here 
cp $EQCHECKPOINT init_0.checkpoint || exit 1

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
	       sleep 10s
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
	       sleep 10s
	    else
	       echo "Continuing $traj-th simulation 20 ps"
	       i-pi $CHECKPOINT &> log_$traj &
	       sleep 10s
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
	       sleep 10s
	    fi
	fi
    # After each traj, we copy the final configuration as the initial configuration
    cp simu_$traj.checkpoint init_$traj.checkpoint || exit 1
done
wait 

