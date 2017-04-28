/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
// #include "petscvec.h"
#include "petsc.h"
#include "petscmat.h"

PetscErrorCode real2complexMat(Mat *A,Mat * r2c);
PetscErrorCode real2complexVec(Vec *in, Vec * out);
