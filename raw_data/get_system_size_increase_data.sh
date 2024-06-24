#!/bin/bash

#SBATCH --job-name=post_processing
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli

python calculate_ph_sp_systemsize_increase.py co2_liquid_phase_systemsize/size_times_*
#python calculate_ph_sp_systemsize_increase.py co2_liquid_phase_systemsize/mode_times_*
