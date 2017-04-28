/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef CONVHULL_H_

#include "petsc.h"
#include <math.h>

int convhull(PetscScalar * ab, PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt * ne, PetscInt offset, PetscInt mu);


#endif