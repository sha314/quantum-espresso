# Ubuntu 24.04

# nvidia-driver-570
# Cuda 12.5
# nvhpc-25-3


# Add NVIDIA repo and install drivers
sudo add-apt-repository ppa:graphics-drivers/ppa -y
sudo apt update
# Install the correct driver (match user-space version)
sudo apt install nvidia-driver-570 -y  # Replace "570" with your target version
sudo reboot

## nvidia-smi shows the maximum CUDA version supported by your driver (e.g., driver 550+ supports CUDA 12.8).
nvidai-smi


## Avoid partial upgrades
sudo apt full-upgrade -y  # Ensures all packages sync


## Pin driver version (if needed):
echo "Package: *nvidia*" | sudo tee /etc/apt/preferences.d/nvidia-pin
echo "Pin: version 570*" | sudo tee -a /etc/apt/preferences.d/nvidia-pin
echo "Pin-Priority: 1001" | sudo tee -a /etc/apt/preferences.d/nvidia-pin
sudo apt update


## Install Cuda toolkit (may not install 12.5)
# sudo apt install nvidia-cuda-toolkit
# nvcc --version



#### https://developer.nvidia.com/cuda-12-5-0-download-archive?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=22.04&target_type=deb_network
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-5
nvcc --version


## Probably nvhpc-25-3-cuda-multi will work?


curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y
sudo apt-get install -y nvhpc-25-3


# curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
# echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
# sudo apt-get update -y
# sudo apt-get install -y nvhpc-25-3-cuda-multi


# Compilers
nvc --version
nvc++ --version
nvfortran --version



# Test the compiler
# Create a file "test.cuf" and add the folloing code there

program hello
  use cudafor
  write(*,*) "CUDA Fortran works with NVHPC 2025!"
end program hello

# Compile and run it
nvfortran -cuda test.cuf -o test
./test



## Environment variables for nvcc, nvc, nvc++, nvfortran


echo 'export NVHPC=/opt/nvidia/hpc_sdk' >> ~/.bashrc
echo 'export PATH=$NVHPC/Linux_x86_64/25.3/compilers/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$NVHPC/Linux_x86_64/25.3/compilers/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export MANPATH=$NVHPC/Linux_x86_64/25.3/compilers/man:$MANPATH' >> ~/.bashrc

source ~/.bashrc


## Environment variables for mpirun

echo 'export HPCX_DIR=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1' >> ~/.bashrc
echo 'export PATH=$HPCX_DIR/ompi/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$HPCX_DIR/ompi/lib:$LD_LIBRARY_PATH' >> ~/.bashrc

source ~/.bashrc


