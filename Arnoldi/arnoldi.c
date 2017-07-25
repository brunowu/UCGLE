/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "arnoldi.h"

/* Compute cyclicly eigenvalue */
PetscErrorCode Arnoldi(com_lsa * com, Mat * A, Vec  *v){

	EPS eps; /* eigensolver context */
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	PetscInt end,first,validated;
	PetscErrorCode ierr;
	/* eigenvalues number is set to 100, can be changed if needed
	   we choosed to fix it because mallocs weren't working properly */
	PetscScalar eigenvalues[500], ei, er;
	PetscReal re,im,vnorm;
	PetscInt eigen_nb,j,i,size, taille;
	uintptr_t one=1;
	Vec initialv,nullv;
	PetscBool flag,data_load,data_export,continuous_export,load_any, aft_flg;
	int exit_type=0;
	int sos_type = 911;
	Vec vecteur_initial;
	PetscViewer viewer;
	PetscInt aft, count=0;

	char  ls_load_path[PETSC_MAX_PATH_LEN];
	PetscBool ls_load, ls_load_any;

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_ls_load",ls_load_path,PETSC_MAX_PATH_LEN,&ls_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_ls_load_any",&ls_load_any);CHKERRQ(ierr);

	if(ls_load&&ls_load_any){
	  ls_load_any=PETSC_FALSE;
	  ls_load=PETSC_TRUE;
	}
	PetscBool need_new_init = PETSC_FALSE, exit = PETSC_FALSE;

	end=0;
	first=1;
	validated=1;

 if(!(ls_load^=ls_load_any)){

	sprintf(load_path,"./arnoldi.bin");
	sprintf(export_path,"./arnoldi.bin");

	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);

	/* create the eigensolver */
	ierr=EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	/* set the matrix operator */
	ierr=EPSSetOperators(eps,*A,PETSC_NULL);CHKERRQ(ierr);
	ierr=EPSSetProblemType(eps, EPS_NHEP);CHKERRQ(ierr);
	/* set options */
	ierr=EPSSetType(eps,EPSARNOLDI);
	ierr=EPSSetFromOptions(eps);CHKERRQ(ierr);
	/* duplicate vector properties */
	ierr=VecDuplicate(*v,&initialv);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&nullv);CHKERRQ(ierr);
	ierr=VecSet(nullv,(PetscScalar)0.0);CHKERRQ(ierr);

	ierr=VecSetRandom(initialv,PETSC_NULL);//initialize initial vector to random
	ierr=VecGetSize(initialv,&size);CHKERRQ(ierr);

	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_eigen",&eigen_nb,&flag);CHKERRQ(ierr);
	if(!flag) eigen_nb=EIGEN_ALL;
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_arnoldi_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_arnoldi_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);

	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_arnoldi_load_any",&load_any);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_arnoldi_cexport",&continuous_export);CHKERRQ(ierr);
	 
	ierr=VecDuplicate(initialv,&vecteur_initial);CHKERRQ(ierr);

	PetscOptionsGetInt(NULL,NULL,"-ArnoldiFT",&aft,&aft_flg);

	if(!aft_flg){
		aft = 20000;
	}
	while(!end && count < aft){
		count ++;
		/*check if the program need to exit */
		if(exit == PETSC_TRUE)
			break;
		/* check if we received an exit message from Father*/
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
		/* check if we received some data from GMRES		*/
		if(!mpi_lsa_com_vec_recv(com, &initialv)){
				VecGetSize(initialv, &taille);
		}

		for(j=0;j<eigen_nb;j++){
			eigenvalues[j]=(PetscScalar)0.0;
		}

		if(data_load&&load_any){
		  load_any=PETSC_FALSE;
		  data_load=PETSC_TRUE;
		}
		ierr = VecAssemblyBegin(initialv);CHKERRQ(ierr);
  	ierr = VecAssemblyEnd(initialv);CHKERRQ(ierr);

		if(!(data_load^=load_any)){
		  ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);
		} else {
		  ierr=readBinaryVecArray(load_path,(int*)one,&initialv);CHKERRQ(ierr);
			data_load=PETSC_FALSE;
			load_any=PETSC_FALSE;
			ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);
		}
		ierr=EPSSolve(eps);CHKERRQ(ierr);
		/* get the number of guessed eigenvalues */
		ierr=EPSGetConverged(eps,&eigen_nb);CHKERRQ(ierr);
		/* send them */
		for(j=0;j<eigen_nb;j++){
			ierr = EPSGetEigenvalue(eps,j,&er,&ei);CHKERRQ(ierr);
			#ifdef PETSC_USE_COMPLEX
			  re=PetscRealPart(er);
			  im=PetscImaginaryPart(er);
			#else
			  re=er;
			  im=ei;
			#endif
			eigenvalues[j]=(PetscScalar)re+PETSC_i*(PetscScalar)im;
		}
		if( eigen_nb != 0){
			mpi_lsa_com_array_send(com, &eigen_nb, eigenvalues);
			if(!mpi_lsa_com_type_recv(com,&exit_type)){
   			if(exit_type==666){
		 			end=1;
		 			mpi_lsa_com_type_send(com,&exit_type);
		 			need_new_init = PETSC_FALSE;
		 			exit = PETSC_TRUE;
					break;
  				}
			}
  	}
	}

if(data_export){
	ierr=writeBinaryVecArray(export_path, 1, &initialv);
}
/* and destroy the eps */
ierr=EPSDestroy(&eps);CHKERRQ(ierr);
ierr=VecDestroy(&initialv);CHKERRQ(ierr);
ierr=VecDestroy(&nullv);CHKERRQ(ierr);
}

else{
	while(!end){
		/*check if the program need to exit */
		if(exit == PETSC_TRUE)
			break;
		/* check if we received an exit message from Father*/
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
   }
	}


return 0;
}
