# Mesoscale Molecular Simulations of Fabry-Pérot Vibrational Strong Coupling

Data for: Mesoscale Molecular Simulations of Fabry-Pérot Vibrational Strong Coupling [arXiv:2403.12282](https://arxiv.org/abs/2403.12282), [JCTC accepted](https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c00349).

The folder **demo/** contains all the necessary files to run a single Fabry-Perot VSC CavMD calculation for cavity #1 in the above publication. One just needs to run the calculation as the normal i-pi calculation with LAMMPS. Here, 36 cavity modes are coupled to 36 molecular grid points, and the molecular forces in different grid points can be evaluated in parallel.

- init.xyz: molecular geometries in 36 grid points + 36 cavity modes.
- input.xml: the section **ffgencavsocket** contains all the information to setup the 1D Fabry-Perot VSC simulation; see also Table 1 in the paper for a one-to-one correspondance of the parameters.


The folder **raw_data/** contains all the input files and raw data to reproduce the figures in the abvoe publication. Please go to **raw_data/plotting/** and run the corresponding python scripts to plot the figures.