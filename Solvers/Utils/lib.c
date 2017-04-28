/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "lib.h"

/* calcul de l'epsilon machine */
PetscReal epsilon(void){
	PetscReal one=(PetscReal)1.0, two=(PetscReal)2.0,temp=(PetscReal)1.0;
	
	while(one+temp>one)
		temp/=two;
		
	return two*temp;
}
