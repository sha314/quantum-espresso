# quantum-espresso
Using Quantum ESPRESSO for DFT calculations. I'll add examples with steps and commands.

# Installation
1. Quantum ESPRESSO (QE)
2. gnuplot
3. CIF to  QE input file converter

# Preparing the input files
1. Pseudo potentials are crucial for calculations. There are two ways, you can have them in "./pseudo" folder. Here "." denotes the directory from where you run the commands. And you'd need to have 'pseudo_dir' defined in &CONTROL

```bat
&CONTROL
    pseudo_dir  = './pseudo'
```

Or 

have the pseudo potential in a different directory and define ESPRESSO_PSEUDO environment variable.

I'm choosing the 2nd way since I don't want to upload the pseudo potentials.


You can download the pseudo potentials from here https://www.materialscloud.org/discover/sssp/table/efficiency . It's a file with ".UPF" extension. For silicon the name is "Si.pbe-n-rrkjus_psl.1.0.0.UPF". I'll not upload the pseudo potentials to github.





2. After running each command, you can check if "JOB DONE" text exists in ".out" file.  



## Commands using OpenMPI and CPU
### for SCF
```bat
mpirun -np 4 pw.x -in si_scf.in > si_scf.out
```

### Non-SCF
```bat
mpirun -np 4 pw.x -in si_nscf.in > si_nscf.out
mpirun -np 4 pw.x -in si_bands.in > si_bands.out
mpirun -np 4 bands.x -in si_bands.pp.in > si_bands.pp.out
mpirun -np 4 dos.x -in si_dos.in > si_dos.out
```

### Warning!!!!!
If you use GPU then you should use "-np 1" or "-n 1"




## Exploring calculated files
GNUplot might be a faster way to make initial plots

```bat
sudo apt install gnuplot
gnuplot plot_bands.gp
```
 

