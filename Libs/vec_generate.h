/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef _READ_MATRIX_H
#define _READ_MATRIX_H

#include "petscmat.h"
#include "petscerror.h"

#endif
PetscErrorCode generate_canoic_vector(PetscInt size, PetscScalar value, PetscInt position, Vec * v);
PetscErrorCode generate_homogene_vector(PetscInt size, PetscInt scal,Vec * v);
PetscErrorCode generate_random_vector(PetscInt size, PetscInt low, PetscInt up,Vec * v);
PetscErrorCode generate_continuzero_vector(PetscInt size, PetscInt start, PetscInt length, PetscScalar value, Vec * v);
PetscErrorCode generate_random_seed_vector(PetscInt size, PetscReal low, PetscReal up, PetscReal seed, Vec * v);
