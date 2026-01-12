# FFLUX_MPI
Most up to date version of the FFLUX code. (Supersedes any other versions).
Fully compatible with DL_POLY's domain decomposition MPI.

# Modules for CSF4
The mpif90 compiler is used. 
```
module load iompi/2020.02
```
After that, run `make hpc`.

# Modules for CSF3
The mpif90 compiler is used.
- If encountered errors in run, delete the -xHost in the Makefile in hpc flags 
```
module load mpi/intel-18.0/openmpi/3.1.4
```

After that, run `make hpc`.


