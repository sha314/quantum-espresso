#!/bin/bash
# Everything is compiled. Only the system compilers are needed to be installed
# Optimized for 7970x Threadripper


##################################################
# Directory for user
##################################################
# Create working directory 
# Change ownership to your user # Give full access to your user
sudo mkdir -p /opt/7970x
sudo chown -R $USER:$USER /opt/7970x
sudo chmod -R 755 /opt/7970x


#################################################
# Set PATH and LIB
#################################################

# better runtime defaults for single-node MPI
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1



export OMP_NUM_THREADS=2

# MPI
export PATH=/opt/7970x/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/opt/7970x/openmpi/lib:$LD_LIBRARY_PATH

# OpenBLAS
export LD_LIBRARY_PATH=/opt/7970x/openblas/lib:$LD_LIBRARY_PATH

# FFTW
export LD_LIBRARY_PATH=/opt/7970x/fftw/lib:$LD_LIBRARY_PATH

# Parallel HDF5
export PATH=/opt/7970x/hdf5/bin:$PATH
export LD_LIBRARY_PATH=/opt/7970x/hdf5/lib:$LD_LIBRARY_PATH

# Quantum ESPRESSO
export PATH=/opt/7970x/qe/bin:$PATH
export LD_LIBRARY_PATH=/opt/7970x/qe/lib:$LD_LIBRARY_PATH



##################################################
# Ryzen 7970X tuning (Zen 4)
##################################################

# 7970X = Zen 4 → use znver4, not znver2
# -Ofast usually performs better for QE than -O3
# keep fPIC for shared libs

export CFLAGS="-Ofast -march=znver4 -mtune=znver4 -fPIC"
export CXXFLAGS="$CFLAGS"
export FCFLAGS="-Ofast -march=znver4 -mtune=znver4 -fPIC -fallow-argument-mismatch"
export FFLAGS="$FCFLAGS"



##################################################
# OpenMPI (same toolchain for everything)
##################################################

    cd /tmp

    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz

    tar xf openmpi-4.1.6.tar.gz

    cd openmpi-4.1.6

    ./configure \
        --prefix=/opt/7970x/openmpi \
        --enable-mpi-fortran \
        CC=gcc \
        CXX=g++ \
        FC=gfortran

    make -j$(nproc)
    make install

    

    export CC=mpicc
    export CXX=mpicxx
    export FC=mpifort


##################################################
# OpenBLAS optimized for EPYC
##################################################

    cd /tmp

    wget https://github.com/xianyi/OpenBLAS/archive/v0.3.26.tar.gz

    tar xf v0.3.26.tar.gz

    cd OpenBLAS-0.3.26

    make -j$(nproc) \
    CC=gcc \
    FC=gfortran \
    TARGET=ZEN \
    USE_OPENMP=1 \
    NO_MPI=1

    make PREFIX=/opt/7970x/openblas install


##################################################
# FFTW3 MPI build
##################################################

    cd /tmp

    wget https://www.fftw.org/fftw-3.3.10.tar.gz

    tar xf fftw-3.3.10.tar.gz

    cd fftw-3.3.10

    ./configure \
        --prefix=/opt/7970x/fftw \
        --enable-mpi \
        --enable-shared \
        CC=mpicc

    make -j$(nproc)

    make install


# ##################################################
# # ScaLAPACK (fixed version)
# ##################################################

#     cd /tmp

#     wget https://github.com/Reference-ScaLAPACK/scalapack/archive/v2.2.1.tar.gz

#     tar xf v2.2.1.tar.gz

#     cd scalapack-2.2.1

#     mkdir build
#     cd build

#     cmake .. \
#     -DCMAKE_INSTALL_PREFIX=/opt/7970x/scalapack \
#     -DCMAKE_C_COMPILER=mpicc \
#     -DCMAKE_Fortran_COMPILER=mpifort \
#     -DBLAS_LIBRARIES=/opt/7970x/openblas/lib/libopenblas.so \
#     -DLAPACK_LIBRARIES=/opt/7970x/openblas/lib/libopenblas.so \
#     -DCMAKE_Fortran_FLAGS="$FCFLAGS"

#     make -j$(nproc)
#     make install



##################################################
# Build Parallel HDF5
##################################################
cd /tmp
git clone https://github.com/HDFGroup/hdf5.git
cd hdf5

git checkout hdf5_2.1.0
mkdir build && cd build

cmake .. \
-DCMAKE_INSTALL_PREFIX=/opt/7970x/hdf5 \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_Fortran_COMPILER=mpifort \
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



##################################################
# Skip scalapack for single CPU
# Quantum ESPRESSO build. 7.5
##################################################
   cd /tmp

   git clone --depth 1 --branch qe-7.5 https://github.com/QEF/q-e.git
   cd q-e

   mkdir build
   cd build

   cmake .. \
   -DCMAKE_INSTALL_PREFIX=/opt/7970x/qe \
   -DCMAKE_C_COMPILER=mpicc \
   -DCMAKE_Fortran_COMPILER=mpifort \
   -DQE_ENABLE_MPI=ON \
   -DBLAS_LIBRARIES=/opt/7970x/openblas/lib/libopenblas.so \
   -DLAPACK_LIBRARIES=/opt/7970x/openblas/lib/libopenblas.so \
   -DFFTW3_ROOT=/opt/7970x/fftw \
   -DHDF5_ROOT=/opt/7970x/hdf5 \
   -DHDF5_INCLUDE_DIRS=/opt/7970x/hdf5/include \
   -DHDF5_LIBRARIES="/opt/7970x/hdf5/lib/libhdf5_fortran.so;/opt/7970x/hdf5/lib/libhdf5.so"

   make -j$(nproc)
   make install



##################################################
# Cleanup
##################################################

    rm -rf /tmp/openmpi*
    rm -rf /tmp/OpenBLAS*
    rm -rf /tmp/fftw*
    rm -rf /tmp/scalapack*
    rm -rf /tmp/hdf5*
    rm -rf /tmp/vscode-server*
    rm -rf /tmp/v*.tar.gz
    rm -rf /tmp/*.tar.gz


