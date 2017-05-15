/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef GMRES_SOLVE_H
#define GMRES_SOLVE_H

#include "petsc.h"
#include "gmres_cycle.h"
#include "gmres_precond.h"
#include "../../Libs/mpi_lsa_com.h"
#include "fgmresimpl.h"
//#include </home/xinzhewu/Petsc/petsc-3.7.0/src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>
//#include <../src/ksp/ksp/impls/gmres/fgmres/fgmresimpl.h>

PetscErrorCode MyKSPSolve_FGMRES(KSP ksp,com_lsa * com);

#endif
