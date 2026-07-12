#!/bin/bash

## CPU-only Phoebe / Quantum ESPRESSO environment
## Compiled with: GNU GCC/gfortran, OpenMPI 4.1.6, UCX 1.21.0, OpenBLAS 0.3.33,
##               ScaLAPACK 2.2.3, FFTW 3.3.10, Parallel HDF5 2.1.0
## Target: Zen 4 (x86_64), MPI+OpenMP, no GPU
##

conda deactivate # deactivate all conda env


export OMP_NUM_THREADS=1

export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_pml=ucx
export UCX_TLS=self,sm,tcp
export UCX_NET_DEVICES=mlx5_0:1


# UCX
export PATH=/opt/ucx/bin:$PATH
export LD_LIBRARY_PATH=/opt/ucx/lib:$LD_LIBRARY_PATH

# MPI
export PATH=/opt/openmpi_gnu/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi_gnu/lib:$LD_LIBRARY_PATH

# OpenBLAS
export LD_LIBRARY_PATH=/opt/openblas/lib64:$LD_LIBRARY_PATH

# FFTW
export LD_LIBRARY_PATH=/opt/fftw/lib:$LD_LIBRARY_PATH

# ScaLAPACK
export LD_LIBRARY_PATH=/opt/scalapack/lib64:$LD_LIBRARY_PATH

# Parallel HDF5
export PATH=/opt/hdf5_2.1.0_gnu/bin:$PATH
export LD_LIBRARY_PATH=/opt/hdf5_2.1.0_gnu/lib:$LD_LIBRARY_PATH


## Phoebe patched QE 7.3 
export QE=/home/$USER/Downloads/PHOEBE/phoebe-quantum-espresso/build_gnu_CPU
export QEX=$QE/build_gnu_CPU/bin
export EPW=$QE/EPW/bin

export PHOEBE=/home/$USER/Downloads/PHOEBE/phoebe
export PHOEBEX=$PHOEBE/build_gnu_CPU


##################################################
### System software 
##################################################


sudo dnf install binutils-devel
sudo dnf install libibverbs-devel rdma-core-devel
sudo dnf install hwloc-devel libevent-devel


sudo chown -R $USER:$USER /opt
sudo chmod -R 755 /opt/

export TEMPORARY=~/Downloads/PHOEBE

##################################################
# UCX 
# To support, InfiniBand, RoCE, or plan to use AMD GPUs for GPU-aware MPI.
##################################################

cd $TEMPORARY

wget https://github.com/openucx/ucx/archive/refs/tags/v1.21.0.tar.gz
tar xf v1.21.0.tar.gz 
cd ucx-1.21.0

./autogen.sh
./configure --prefix=/opt/ucx \
    --enable-mt \
    --without-go \
    --disable-tests \
    CC=gcc CXX=g++

make -j$(nproc)
make install





##################################################
# OpenMPI
##################################################
cd $TEMPORARY
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz

tar xf openmpi-4.1.6.tar.gz

cd openmpi-4.1.6

./configure --prefix=/opt/openmpi_gnu \
    --with-ucx=/opt/ucx \
    --without-ucc \
    --enable-mpi-thread-multiple \
    --enable-mpi-cxx \
    --enable-shared \
    CC=gcc CXX=g++ FC=gfortran

make -j$(nproc)
make install



##################################################
# OpenBLAS 
##################################################

cd $TEMPORARY

wget https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v0.3.33.tar.gz

tar xf v0.3.33.tar.gz 
cd OpenBLAS-0.3.33/

mkdir build && cd build

cmake .. \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_FLAGS="-O3 -march=znver4 -mtune=znver4 -fPIC" \
    -DCMAKE_Fortran_FLAGS="-O3 -march=znver4 -mtune=znver4 -fPIC" \
    -DTARGET=ZEN4 \
    -DUSE_OPENMP=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=/opt/openblas


make -j$(nproc)
make install




##################################################
# FFTW3 MPI build
##################################################

cd $TEMPORARY
wget https://www.fftw.org/fftw-3.3.10.tar.gz

tar xf fftw-3.3.10.tar.gz

cd fftw-3.3.10

./configure \
    --prefix=/opt/fftw \
    --enable-mpi \
    --enable-shared \
    --enable-threads \
    --enable-openmp \
    CC=mpicc \
    MPICC=mpicc \
    FC=mpif90 \
    CFLAGS="-O3 -fPIC"

make -j$(nproc)

make install




##################################################
# ScaLAPACK 
##################################################

cd $TEMPORARY

wget https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v2.2.3.tar.gz

tar xf v2.2.3.tar.gz 

cd scalapack-2.2.3

mkdir build
cd build

cmake .. \
    -DCMAKE_INSTALL_PREFIX=/opt/scalapack \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_C_FLAGS="-O3 -fPIC" \
    -DCMAKE_Fortran_FLAGS="-O3 -fPIC" \
    -DBLAS_LIBRARIES=/opt/openblas/lib64/libopenblas.so \
    -DLAPACK_LIBRARIES=/opt/openblas/lib64/libopenblas.so 

make -j$(nproc)
make install



##################################################
# Build Parallel HDF5
##################################################
cd $TEMPORARY
git clone https://github.com/HDFGroup/hdf5.git
cd hdf5

git checkout hdf5_2.1.0
mkdir build && cd build



cmake .. \
    -DCMAKE_INSTALL_PREFIX=/opt/hdf5_2.1.0_gnu \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_C_FLAGS="-O3 -fPIC" \
    -DCMAKE_Fortran_FLAGS="-O3 -fPIC" \
    -DHDF5_BUILD_FORTRAN=ON \
    -DHDF5_ENABLE_PARALLEL=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_TOOLS=ON \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF


make -j$(nproc)
make install

# Check if you have the .mod files for the Fortran modules needed to build QE
find /opt/hdf5_2.1.0_gnu  -name "*.mod"


##################################################
# Phoebe Patched Quantum Espresso 7.3
##################################################
cd $TEMPORARY
git clone https://github.com/phoebe-team/phoebe-quantum-espresso.git
cd phoebe-quantum-espresso
git checkout phoebe-qe-7.3
git submodule update --init --recursive



mkdir build_gnu_CPU
cd build_gnu_CPU

cmake .. \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DCMAKE_BUILD_TYPE=Release \
    -DQE_ENABLE_MPI=ON \
    -DQE_ENABLE_SCALAPACK=ON \
    -DQE_ENABLE_HDF5=ON  \
    -DQE_FFTW_VENDOR=Internal \
    -DHDF5_ROOT=/opt/hdf5_2.1.0_gnu

make -j$(nproc)




##################################################
# Phoebe v1.1.0
##################################################
cd $TEMPORARY
git clone https://github.com/phoebe-team/phoebe.git
cd phoebe
git checkout v1.1.0
git submodule update --init --recursive


mkdir build_gnu_CPU
cd build_gnu_CPU

cmake .. \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_STANDARD=20 \
    -DMPI_AVAIL=ON \
    -DOMP_AVAIL=ON \
    -DKokkos_ENABLE_OPENMP=ON \
    -DHDF5_AVAIL=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS="-lgfortran -lquadmath"


make -j$(nproc)





