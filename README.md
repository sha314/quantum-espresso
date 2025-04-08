# quantum-espresso
Using Quantum ESPRESSO for DFT calculations. I'll add examples with steps and commands.

# Installation
1. Quantum ESPRESSO (QE)
2. gnuplot
3. CIF to  QE input file converter

# Preparing
1. Make sure you have the pseudo potentials in "./pseudo" folder. Here "." denotes the directory from where you run the commands. You can download the pseudo potentials from here https://www.materialscloud.org/discover/sssp/table/efficiency . It's a file with ".UPF" extension. For silicon the name is "Si.pbe-n-rrkjus_psl.1.0.0.UPF". I'll not upload the pseudo potentials to github.


2. After running each command, you can check if "JOB DONE" text exists in ".out" file.  



# Commands using OpenMPI
### for SCF
mpirun -np 4 pw.x -in Si/si_scf.in > si_scf.out

### Non-SCF
mpirun -np 4 pw.x -in Si/si_nscf.in > si_nscf.out

mpirun -np 4 pw.x -in Si/si_bands.in > si_bands.out

mpirun -np 4 bands.x -inp Si/si_bands.pp.in > si_bands.out

mpirun -np 4 dos.x -in Si/si_dos.in > si_dos.out






