# cavmd_examples_fpvsc

Please stay tuned, the data will be uploaded upon the publication of the paper: Mesoscale Molecular Simulations of Fabry-Pérot Vibrational Strong Coupling [arXiv:2403.12282, JCTC accepted](https://arxiv.org/abs/2403.12282).

Currently, the folder **demo/** contains all the necessary files to run a single Fabry-Perot VSC CavMD calculation for cavity #1 in the above publication. One just needs to run the calculation as the normal i-pi calculation with LAMMPS. Here, 36 cavity modes are coupled to 36 molecular grid points, and the molecular forces in different grid points can be evaluated in parallel.

init.xyz: molecular geometries in 36 grid points + 36 cavity modes.
input.xml: the section **ffgencavsocket** contains all the information to setup the 1D Fabry-Perot VSC simulation; see also Table 1 in the paper for a one-to-one correspondance of the parameters.

All the other data will be organized and uploaded once the paper is online in JCTC.
