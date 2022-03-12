# CHORUS-MHD
This code uses the high-order Spectral Difference method with Unstaggered Constrained Transport to perform compressible magnetohydrodynamic simulations.
This method can preserve the divergence-free constraint exactly and globally.
For details about the numerical methods, please refer to the journal publication
https://doi.org/10.1080/10618562.2022.2042272

Systems that have been tested on:
Linux (various flavors/distros), 64 bit (x86), with a fortran compiler (gfortran or ifort) and a MPI compiler (openmpi or mpich)

Commands to run this code:
make

cd rundir

mpirun -np number_of_CPU ./sd2d_mhd

Commands to clean .o and executable files:
make clean

Parallel run:
please install METIS library for mesh partitioning
The most updated METIS library can be downloaded from
http://glaros.dtc.umn.edu/gkhome/metis/metis/download
The file for mesh partitioning an .mts file
Use the metis command
mpmetis filename.mts number_of_CPU
to partition the mesh for parallel computing

For further questions, please get contact with Kuangxu (Scott) Chen at Clarkson University: kuchen@clarkson.edu
