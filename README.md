# metadynamics-tools

xparamber1:
AMBER parametrization automation
Script that automates the process of parameterizing small ligand molecules in the AMBER forcefield using the AmberTools package.

walkermaker:
Tools for generating walkers
These tools help with generating multiple starting positions for a ligand to be used in a metadynamics simulation.


restminimizer:
Scripts for running em with restrains
The ndxmaker.sh program generates an .ndx file, which is then used to create a .posre file.
The runrestrained_em.sh program performs energy minimization for multiple walkers using a .posre file.
