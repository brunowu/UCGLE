# UCGLE method (Unite and Conquer GMRES/LS-ERAM method)


This manuel comprises the source code, datasets, and build instructions that can be used to reproduce the results of our UCGLE method and GMRES_LS method, we explain how to compile, install and run our UGCLE method codes, all the results of this paper can be reproduced with the help of this description.

# Description

## Check-list


• Algorithm: UCGLE method and GMRES_LS method

• Program: C MPI codes

• Compilation: C compiler (gcc version 4.9.1 is tested) 

• Binary: MPI executable

• Data set: Public available matrix file from Matrix Market and generated matrix files

• Run-time environment: Linux OS (Redhat 6.3 is tested) and Mac OS X (higer than Mac OS X 10.9 are tested)

• Hardware: Any Intel CPU (Intel(R) Xeon(R) CPU E5-2650 v2 @ 2.60GHz and Intel Core i7 @ 2.70GHZ are tested)

• Output: Execution time to solution and interation number

• Experiment workflow: Install PETSc and SLEPc, clone UCGLE codes, compile the codes, run the binary, observe the results

• Experiment customization: numberofGMRESMPIprocessors, number of Arnoldi MPI processors, standard parameters of GMRES
and Arnoldi, Least Square polynomial degree· · ·

• Publicly available?: Yes

## Software Dependencies

The UCGLE method evaluation requires the Intel C compiler, such as gcc, icpc, clang etc. Our codes are based on the scien c libraries PETSc and SLEPc. The newest stable version of PETSC is 3.7.5, of SLEPc is 3.7, but we developped our software with PETSc 3.6.4 and SLEPc 3.6, so if you want run this application with higher versions of libraries, it is necessary to adjust the codes with the new changes, espacially the group of PetscOptions functions. BLAS and LAPACK are need but can be automatically downloaded by PETSc. It is also necessary to support MPI (OpenMPI or MPICH), with can also be automatically installed by PETSc with the con gure  ags. It is very important to well con gure the installation of PETSc (make sure to turn off the debugging mode of PETSc with the  ag –with-debugging=0):

## Datasets

The Datasets include mainly the tested matrices.

• The industrial matrices downloaded from the site Matrix Market, some of them converge with our method, some
  don’t. You can use the example matrix utm300 kept in the directory UCGLE/data.
  
• The generated matrices matLine and matBlock presented

The method is tested on Redhat 6.4 and Mac OS Sierra, and is also expected to run correctly under other Linux based systems.

# Installation

The installation steps:

• Clone UCGLE code to the local machine (the repository will be provided if this paper is accepted)

• Use GNU make to generate the executable binary hyperh

It will come soon the manuel...
