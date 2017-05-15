/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef TRIEP_H_

#include "petsc.h"
#include <math.h>
#include <stdlib.h>

PetscErrorCode tri(PetscScalar * vp, PetscInt nValues, PetscInt * ch_signe);

PetscErrorCode epurer(PetscScalar * vp, PetscInt * nValues);

void swap(PetscScalar *a, PetscScalar *b);

void sort(PetscScalar arr[], int beg, int end);

#endif
