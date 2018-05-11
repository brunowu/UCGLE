/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "read_matrix.h"
#include <interface/C/c_wrapper.h>
PetscErrorCode read_matrix_vector(Mat * A, Vec * v, int *comm){
	char filea[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flaga,flagb, flgsmg2s;
	PetscViewer fd;
	PetscInt size,sizea;
	PetscScalar scal;
	PetscReal vnorm;

	ierr= PetscOptionsHasName(NULL,NULL,"-smg2s",&flgsmg2s);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);

	PetscLogDouble st, ed;

	PetscTime(&st);
	
	if (!flaga) {
		if(!flgsmg2s){
		        sprintf(err,"Error : mfile is not properly set -> %s\n",filea);
	                SETERRQ(*comm,(PetscErrorCode) 83,err);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD,"Remainder ]> USING SMG2S to generate test matrix\n");
			mat_generate(A, *comm);
		}
	}
	else{
		/* read matrix file */
	        PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
       		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
        	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
        	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);

	}

	PetscTime(&ed);

	PetscPrintf(*comm, "In group: The test matrix generation/loading total time is %f \n",ed-st);
	/* read matrix file */
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);

	if (!flagb) {
		/* the user did not provide a vector, so generate it*/
		generate_random_seed_vector(size, -10.0, 10.0, 0, v);
		VecNorm(*v, NORM_2, &vnorm);
		VecScale(*v, 0.01/vnorm);		
		PetscPrintf(PETSC_COMM_WORLD,"Generated right hand side matrix b\n");		
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",fileb);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileb,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		VecCreate(PETSC_COMM_WORLD,v);
		ierr=VecLoad(*v,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

	return 0;
}
