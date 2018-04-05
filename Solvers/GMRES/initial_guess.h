/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef GMRES_PRECOND_H
#define GMRES_PRECOND_H

#include "fgmresimpl.h"
// #include "../Arnoldi/arnoldi.h"
//#include <../src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>
//#include </home/xinzhewu/Petsc/petsc-3.7.0/src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>
#include "petsc.h"
#include "../Utils/lib.h"
#include "../../Libs/mpi_lsa_com.h"
#include <unistd.h>

#ifndef EIGEN_ALL
#define EIGEN_ALL 100
#endif

#ifndef LS_POWER
#define LS_POWER 15
#endif

#ifndef LS_LATENCY
#define LS_LATENCY 2
#endif

PetscErrorCode LSInitialGuess(com_lsa * com, KSP ksp);
//PetscErrorCode InitVecPrecond(com_lsa * com, Mat Amat, Vec b, Vec x_unprecond, Vec * x_initial);
#endif
