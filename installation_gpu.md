
# This document is a guide to compile Quantum ESPRESSO for GPU (RTX 4090) on Ubuntu 24.04

## Step 1: GPU Drivers (570), CUDA Toolkit (12.5), NVIDIA HPC SDK (25.3)

Run the commands in 'nvidia-drivers.sh' file


## Step 2: Other packages (FFT, cmake, git)

For GPU only built,
You do NOT need libfftw3 or libfftw3-mpi-dev if:
You're compiling with GPU support and using CUFFT (NVIDIA's FFT library).
You're using cmake with -DQE_ENABLE_CUDA=ON and it detects CUFFT.

But for hybrid built (CPU+GPU), you do need it.

```bat
sudo apt install git cmake
sudo apt-get install libfftw3-dev libfftw3-mpi-dev
```


## Step 3: Quantum ESPRESSO (7.4) with GPU support
Start form home directory

```bat
cd ~
wget https://www.quantum-espresso.org/rdm-download/488/v7-4-1/4bb7e3409081909b8600a33c332a02e5/qe-7.4.1-ReleasePack.tar.gz
tar xvf qe-7.4.1-ReleasePack.tar.gz 
cd qe-7.4.1
mkdir build
cd build
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DQE_ENABLE_CUDA=ON \
  -DBLAS_LIBRARIES="-lblas" \
  -DLAPACK_LIBRARIES="-llapack"
make pw -j 6
make pw ph pp neb pwcond -j 6
ls ./bin/
```

The make command will build executables in "build/bin/" folder.



## Step 3B: Quantum ESPRESSO (7.4) with Hybrid ()
todo

## Setp 4: Setting up


```bat
/opt/qe741/
├── bin/
│   ├── pw.x
│   ├── ph.x
│   └── ...
├── lib/             ← optional (if used)
├── include/         ← optional .mod files (if needed for linking)
├── pseudo/          ← if you use UPF pseudopotentials
└── examples/        ← optional
```

I'll have only bin, lib and include folder. I'll store pseudo potential in different location


```bat
sudo mkdir -p /opt/qe741/bin
sudo mkdir -p /opt/qe741/lib
sudo mkdir -p /opt/qe741/include
```

Copy the files

```bat
cd ~/qe-7.4.1
sudo cp -r build/bin/ /opt/qe741/
sudo cp -r build/lib/ /opt/qe741/
```

Find the .mod files and put them in include 

```bat
find build/ -name "*.mod" -exec sudo cp {} /opt/qe741/include/ \;
```

Change permission and make them executable
```bat
sudo chown -R root:root /opt/qe741/
sudo chmod -R a+rx /opt/qe741/
```


### Environment variables


```bat
echo 'export PATH=/opt/qe741/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/opt/qe741/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export CPATH=/opt/qe741/include:$CPATH' >> ~/.bashrc
source ~/.bashrc
```


And ~/pseudo will be the Pseudo potential directory (.UPF files).

```bat
echo 'export ESPRESSO_PSEUDO=~/pseudo' >> ~/.bashrc
source ~/.bashrc
```
No need to use full paths in the input files — QE will search $ESPRESSO_PSEUDO.




## Step 5. Testing

Use a GPU-aware executable (pw.x must be GPU-enabled).
```bat
ldd ./pw.x | grep cuda
```
You should see CUDA libraries like libcufft, libcublas, or libcuda






