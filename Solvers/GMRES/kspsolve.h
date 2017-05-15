/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef MY_KSP_SOLVE_H
#define MY_KSP_SOLVE_H
#include "petsc.h"
#include "../../Libs/mpi_lsa_com.h"
#include "gmres_solve.h"
#include <petsc/private/kspimpl.h>

PetscErrorCode MyKSPSolve(KSP ksp,Vec b,Vec x,com_lsa * com);


#endif
