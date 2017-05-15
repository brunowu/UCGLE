/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "petsc.h"
#include "petscvec.h"

PetscErrorCode writeBinaryScalarArray(const char * name, int nb, PetscScalar * array);
PetscErrorCode writeBinaryVecArray(const char * name, int nb, Vec * array);

PetscErrorCode readBinaryScalarArray(const char * name, int * nb, PetscScalar * array);
PetscErrorCode readBinaryVecArray(const char * name, int * nb, Vec * array);

PetscErrorCode getFileSize(const char * name, long * size);
