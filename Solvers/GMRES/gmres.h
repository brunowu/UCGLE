/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef GMRES_H
#define GMRES_H

//#include "fgmresimpl.h"
#include "petsc.h"
#include "petscvec.h"
#include "petscksp.h"
#include "mpi.h"
#include "../../Libs/mpi_lsa_com.h"
#include "gmres_solve.h"
#include "kspsolve.h"
#include <time.h>
#include <unistd.h>
//#include "gmres_precond.c"

#include "../../Libs/vec_generate.h"

#include "gmres_precond.h"
#include "initial_guess.h"


PetscErrorCode launchGMRES(com_lsa * com, Vec * b, Mat * A);

#endif
