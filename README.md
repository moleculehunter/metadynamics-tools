# metadynamics-tools


**xparamber1**:

*AMBER parametrization automation*

Script that automates the process of parameterizing small ligand molecules in th
e AMBER forcefield using the AmberTools package.

---

**walkermaker**:

*Tools for generating walkers*

These tools help with generating multiple starting positions for a ligand to be
used in a metadynamics simulation.

---

**restminimizer**:

*Scripts for running em with restrains*

The ndxmaker.sh program generates an .ndx file, which is then used to create a .
posre file.

The runrestrained_em.sh program performs energy minimization for multiple walker
s using a .posre file.

---

**reboxing**

*Move your walkers into new box*

fixpbc.sh - repair Ligand destroyed by PBC

solver.sh - Change the type and size of box, and then solvate system

ionizer.sh - Neutralize system with K and CL ions
