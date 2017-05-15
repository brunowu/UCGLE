/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef GMRES_CYCLE_H
#define GMRES_CYCLE_H

/* just one include necessary */
#include "fgmresimpl.h"
//#include <../src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>
//#include </home/xinzhewu/Petsc/petsc-3.7.0/src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>
#include "petsc.h"
#include "petsc/private/kspimpl.h"


PetscErrorCode MyKSPFGMRESCycle(PetscInt *itcount,KSP ksp);
PetscErrorCode MyKSPFGMRESGetNewVectors(KSP,PetscInt);
PetscErrorCode MyKSPFGMRESUpdateHessenberg(KSP,PetscInt,PetscBool,PetscReal *);
PetscErrorCode MyKSPFGMRESBuildSoln(PetscScalar*,Vec,Vec,KSP,PetscInt);

#endif
