#!/bin/bash

#SBATCH --job-name=post_processing
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli

# photon dynamics
#python calculate_ph_dynamics.py co2_liquid_phase/E0_5e-5_excited/ 36
# dipole moment spectrum out cavity
#python calculate_dipole_traj.py co2_liquid_phase/pi_E0_0e-4/ 
# photon spectrum
#python calculate_ph_sp.py co2_liquid_phase/pi_E0_5e-5/ 36
#python calculate_ph_sp.py co2_liquid_phase_ud/E0_5e-5/ 36
python calculate_ph_sp.py co2_liquid_phase_ud/cavlifetime_300/ 36
#python calculate_ph_sp.py co2_liquid_phase/E0_5e-5_point/ 36
#python calculate_ph_sp.py co2_liquid_phase/E0_5e-5_2d/ 100
#python calculate_ph_sp.py co2_liquid_phase/E0_5e-5_108modes/ 108
# photon spectrum for a single mode
#python calculate_ph_sp_singlemode.py co2_liquid_phase/E0_5e-5_singlemode/
#python calculate_ph_sp_singlemode.py co2_liquid_phase/pi_E0_5e-5_singlemode/
#python calculate_ph_sp_singlemode.py co2_liquid_phase/E0_5e-5_point_singlemode/
