PWTK is written in plain Tcl. Hence, no compilation is needed. It it's required to have the following installed: (for Ubuntu systems)

### Installation

```
    sudo apt update
    sudo apt install tcl tcllib tcl-tclreadline tcl-thread

    cd ~
    wget "https://pwtk.ijs.si/download/pwtk-3.2.tar.gz"
    cd /opt
    sudo tar -xvf ~/pwtk-3.2.tar.gz
```

### Setting environment
```
    echo "export PATH=/opt/pwtk-3.2:\$PATH" >> ~/.bashrc
    source .bashrc
```

### QE to PWTK
The following command will create a pwtk script for scf/bands computation. good for getting started

```
qe2pwtk -p pw.x scf.in > scf.tcl
qe2pwtk -p pw.x bands.in > bands.tcl
qe2pwtk -p bands.x bands.pp.in > bands.pp.tcl
``` 



### Configuration
To run a PWTK script with MPI or GPU configuration
```
    pwtk scf.tcl prefix="mpirun -np 2" postfix="-npool 6"
```
Or you can edit `~/.pwtk/pwtk.tcl` file and include run configuration

```
# ~/.pwtk/pwtk.tcl
bin_dir   /opt/qe741/bin      ;# path to GPU pw.x
prefix    mpirun -np 2        ;# default MPI launcher
postfix   -npool 2
```

It will run the `pw.x` in the following format
```
prefix pw.x postfix -in si.scf.in > si.scf.out
```




