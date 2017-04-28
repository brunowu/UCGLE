/*
 This file is part of software for the implementation of UCGLE method, under the supervision of Serge G. Petiton
 <serge.petiton@univ-lille1.fr>.
 
 Copyright (C) 2011—. Pierre-Yves AQUILANTI and Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr> in Maison de la Simulation. 
 All rights reserved.
 
 Permission to use, copy, modify and distribute this software for personal and educational use is hereby granted
 without fee, provided that the above copyright notice appears in all copies and that both that copyright notice 
 and this permission notice appear in supporting documentation, and that the names of all authors are not used in 
 advertising or publicity pertaining to distribution of the software without specific, written prior permission. 
 Xinzhe WU and the author make no representations about the suitability of this software for any purpose. It is 
 provided "as is" without express or implied warranty.
 
 You should have received a copy of the GNU Lesser General Public License along with UCGLE.  If not, see 
 <http://www.gnu.org/licenses/>.

 For more information, contact with Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr>.
 
 */

 #include "read_matrix.h"

PetscReal uniform_distribution(PetscReal rangeLow, PetscReal rangeHigh) {
  PetscReal myRand = rand()/(1.0 + RAND_MAX);
  PetscReal range = rangeHigh - rangeLow;
  PetscReal myRand_scaled = (myRand * range) + rangeLow;
  return myRand_scaled;
}

PetscErrorCode generate_canoic_vector(PetscInt size, PetscScalar value, PetscInt position, Vec * v){
	PetscErrorCode ierr;

	VecCreate(PETSC_COMM_WORLD,v);		// we create a PetscObject Vec
	VecSetSizes(*v,PETSC_DECIDE,size); 	// we set it size to size according to size of the matrix A
	VecSetFromOptions(*v);				// setting the object Vec
	VecSetValue(*v,position,value,INSERT_VALUES);
	VecAssemblyBegin(*v);
	VecAssemblyEnd(*v);
	return ierr;
}
PetscErrorCode generate_homogene_vector(PetscInt size, PetscInt scal,Vec * v){

	PetscErrorCode ierr;

	VecCreate(PETSC_COMM_WORLD,v);		// we create a PetscObject Vec
	VecSetSizes(*v,PETSC_DECIDE,size); 	// we set it size to size according to size of the matrix A
	VecSetFromOptions(*v);				// setting the object Vec
	VecSet(*v,scal);
  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
	return ierr;
}

PetscErrorCode generate_random_vector(PetscInt size, PetscInt low, PetscInt up,Vec * v){
	PetscErrorCode ierr;
	PetscRandom    rnd;

	PetscRandomCreate(PETSC_COMM_SELF,&rnd);
	PetscRandomSetType(rnd,PETSCRAND48);
	PetscRandomSetInterval(rnd, low, up);

	VecCreate(PETSC_COMM_WORLD,v);		// we create a PetscObject Vec
	VecSetSizes(*v,PETSC_DECIDE,size); 	// we set it size to size according to size of the matrix A
	VecSetFromOptions(*v);				// setting the object Vec
	VecSetRandom(*v,rnd); 					// we set all the elements of the vector to be random
	PetscRandomDestroy(&rnd);
  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
	return ierr;
}

PetscErrorCode generate_continuzero_vector(PetscInt size, PetscInt start, PetscInt length, PetscScalar value, Vec * v){
	PetscErrorCode ierr;
	PetscInt *row, i;
	PetscScalar *val;

	PetscMalloc1(length, &row);
	PetscMalloc1(length, &val);

	for(i = 0; i < length; i++){
		row[i] = start + i;
		val[i] = value;
	}

	VecCreate(PETSC_COMM_WORLD,v);		// we create a PetscObject Vec
	VecSetSizes(*v,PETSC_DECIDE,size); 	// we set it size to size according to size of the matrix A
	VecSetFromOptions(*v);				// setting the object Vec
  VecSetValues(*v,length, row,val,INSERT_VALUES);
  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
  return ierr;
}

PetscErrorCode generate_random_seed_vector(PetscInt size, PetscReal low, PetscReal up, PetscReal seed, Vec * v){
	PetscErrorCode ierr;
	PetscInt   i;
	PetscScalar  tmp;

	VecCreate(PETSC_COMM_WORLD, v);
	VecSetSizes(*v, PETSC_DECIDE, size);
	VecSetFromOptions(*v);

	srand(seed);

	for(i = 0; i < size; i++){
		tmp = uniform_distribution(low,up);
		VecSetValues(*v, 1, &i, &tmp, INSERT_VALUES);
	}

  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
  return ierr;
}
