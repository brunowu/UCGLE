/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef LSQR_H
#define LSQR_H 

#include "petsc.h"
#include "../../Libs/mpi_lsa.h"
#include "../../Libs/mpi_lsa_com.h"
#include "../../Libs/data_rw.h"
#include "../Utils/triep.h"
#include "../Utils/lib.h"
#include "../Utils/convhul.h"
#include "../Utils/ellipse.h"
#include "./precond.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef EIGEN_MIN
#define EIGEN_MIN 5
#endif

#ifndef EIGEN_ALL
#define EIGEN_ALL 20	
#endif

#ifndef EIGEN_MAX
#define EIGEN_MAX 20	
#endif

PetscErrorCode LSQR(com_lsa * com, int * vector_size);

#endif
