/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef ARNOLDI_H
#define ARNOLDI_H

#include "slepceps.h"
#include <unistd.h>
#include "../../Libs/mpi_lsa.h"
#include "../../Libs/mpi_lsa_com.h"
#include "../../Libs/data_rw.h"

#ifndef EIGEN_ALL
#define EIGEN_ALL 20
#endif

PetscErrorCode Arnoldi(com_lsa * com, Mat * A, Vec * v);

#endif
