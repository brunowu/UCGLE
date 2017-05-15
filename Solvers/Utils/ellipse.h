/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef ELLIPSE_H_

#include "petsc.h"
#include <math.h>
#include <stdlib.h>

PetscErrorCode ellipse(PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt mu, PetscReal * co, PetscReal * ao2, PetscReal * do2, PetscReal * dr, PetscInt * info);

#endif

