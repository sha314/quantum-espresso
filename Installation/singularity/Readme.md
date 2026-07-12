# Summary



# Commands

singularity build --fakeroot 

# NIX containers


/mnt/nfsdata1/images/qe-mpi-epyc.sif 
/mnt/nfsdata1/images/qe-mpi-epyc-062526.sif 


/mnt/nfsdata1/images/mxr2540/qe731mpich.sif
/mnt/nfsdata1/images/mxr2540/qe-mpi-epyc-v10.sif  
/mnt/nfsdata1/images/mxr2540/qe-mpi-epyc-v3.2.sif  
/mnt/nfsdata1/images/mxr2540/qeplotter.sif  
/mnt/nfsdata1/images/mxr2540/testMPI2.sif  
/mnt/nfsdata1/images/mxr2540/testMPI.sif


# Jul 11, 2026
# With UCX support for infiniband
singularity build --fakeroot /mnt/nfsdata1/images/mxr2540/phoebe-qe-mpi-epyc-v11.sif phoebe-qe-mpi-epyc-v11.def  > log1 &
singularity build --fakeroot /mnt/nfsdata1/images/mxr2540/qe-mpi-epyc-v11.sif qe-mpi-epyc-v11.def  > log2 &




