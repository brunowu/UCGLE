/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "main.h"

static const char help[] = "Solve";

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char ** argv){
	com_lsa com;
	int result_code,vsize;
	Mat A,Af; 	
	Vec v,vf;
	PetscErrorCode ierr;
	PetscInt nols;
	PetscBool flag;
	char  tmp_path[PETSC_MAX_PATH_LEN];
	char cCurrentPath[FILENAME_MAX];
	char path[FILENAME_MAX];
	
	time_t firstTime, secondTime;	
	mode_t process_mask = umask(0);
	/* init of MPI and MPI Utils */
	MPI_Init(&argc,&argv);

///	com.com_world=MPI_COMM_WORLD;

	mpi_lsa_init(argc,argv,&com);

///	com.com_world=MPI_COMM_WORLD;

	MPI_Barrier(MPI_COMM_WORLD);

	/* Launch PETSc */
	PETSC_COMM_WORLD=com.com_group; // give group communicator as main petsc communicator
	ierr=SlepcInitialize(&argc,&argv,(char *)0,help);	CHKERRQ(ierr);

//	PetscPrintf(com.com_world,"]> Initializing PETSc/SLEPc\n");

/*	VecCreate()*/

	/* programm options */
	ierr=PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the GMRES/LS-Arnoldi solver", "KSP");CHKERRQ(ierr);

	/* create temp directories */

//	PetscPrintf(com.com_world,"]> Checking temporary directory structure\n");

	/* See if the user choosed a tmp directory in order to unpack the matrix file for each node group */
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-tmp_path",tmp_path,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
	if (!flag) {
		sprintf(tmp_path,"matrices_tmp");
	}


	/* Get the directory name */
 	if (!getcwd(cCurrentPath, sizeof(cCurrentPath))){
			SETERRQ(PETSC_COMM_WORLD,(PetscErrorCode) 83,"cannot write\n");
  }
	/* Create the directory for temp files */
	sprintf(path,"%s/%s",cCurrentPath,tmp_path);
	result_code = mkdir(path, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);


	/* Now create a temp directory for each group */
	MPI_Barrier(MPI_COMM_WORLD);

//	PetscPrintf(com.com_world,"]> Creating directory : %s\n]> Creating the directory structure for each group...\n",path);

	sprintf(tmp_path,"%s",path);
	MPI_Barrier(MPI_COMM_WORLD);


	/* get the names */
	if(com.color_group==GMRES_GR){
		sprintf(path,"%s/gmres",tmp_path);
	}else if(com.color_group==FATHER_GR){
		sprintf(path,"%s/father",tmp_path);
	}else if(com.color_group==ARNOLDI_GR){
		sprintf(path,"%s/arnoldi",tmp_path);
	}else if(com.color_group==LS_GR){
		sprintf(path,"%s/ls",tmp_path);
	}

	/* set permissions */

//	PetscPrintf(com.com_group,"]> Creating directory : %s\n",path);

	result_code = mkdir(path, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
	setenv("PETSC_TMP",path,1);

	/* now read data files */

	MPI_Barrier(MPI_COMM_WORLD);

//	PetscPrintf(com.com_world,"]> Opening input data files\n");

	MPI_Barrier(MPI_COMM_WORLD);
	
	//ierr = MatrCeate(MPI_Comm*) &com.com_world, &A);CHKERRQ(ierr);	// creating Mat A before loading it from the file 
	//ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	read_matrix_vector(&A, &v,(MPI_Comm*) &com.com_world);
	Af=A;
	vf=v;


	PetscOptionsName("-Aclx","Convert complex matrix to real","PetscAclx",&flag);
	if (flag) {
		ierr= real2complexMat(&A,&Af);CHKERRQ(ierr);
		ierr= real2complexVec(&v, &vf);CHKERRQ(ierr);
	}

	A=Af;
	vf=v;


	ierr=PetscOptionsEnd();CHKERRQ(ierr);

	/* now that we got the data, groups will exchange information on data repartition */

//	PetscPrintf(com.com_world,"]> Now exchanging vector sizes...\n");

	mpi_lsa_com_vecsize_init(&com,&v);

	MPI_Barrier(MPI_COMM_WORLD);

///	PetscPrintf(com.com_world,"]> Vector Exchange OK   ...\n");

	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_nopc",&nols,&flag);CHKERRQ(ierr);
	if(!flag) nols=1;
	else if(nols==0) PetscPrintf(com.com_world,"]> not using LSA preconditionning %d\n",nols);
 	time(&firstTime);

	/*create viewers for each group if asked for*/
	PetscOptionsHasName(NULL,PETSC_NULL,"-output_file",&flag);
	if(flag){
	      if(com.color_group==GMRES_GR){
		      //PetscPrintf(com.com_group,"]> Launching GMRES\n");
		      launchGMRES(&com, &v, & A);
	      }else if(com.color_group==ARNOLDI_GR && nols!=0){
		      //PetscPrintf(com.com_group,"]> Launching Arnoldi\n");
		      Arnoldi(&com,&A,&v)	;
	      }	else if(com.color_group==FATHER_GR && nols!=0){
		      //PetscPrintf(com.com_group,"]> Launching Father\n");
		      Father(&com,&v);
	      } else if(com.color_group==LS_GR && nols!=0){
		      VecGetSize(v,&vsize);
		      //PetscPrintf(com.com_group,"]> Launching LS\n");
		      LSQR(&com,&vsize);
	      }
	}



	if(com.color_group==GMRES_GR){
		PetscPrintf(com.com_group,"]> Launching GMRES\n");
		launchGMRES(&com, &v, &A);
	}else if(com.color_group==ARNOLDI_GR && nols!=0){
		PetscPrintf(com.com_group,"]> Launching Arnoldi\n");
		Arnoldi(&com,&A,&v)	;
	}	else if(com.color_group==FATHER_GR && nols!=0){
		PetscPrintf(com.com_group,"]> Launching Father\n");
		Father(&com,&v);
	} else if(com.color_group==LS_GR && nols!=0){
		VecGetSize(v,&vsize);
		PetscPrintf(com.com_group,"]> Launching LS\n");
		LSQR(&com,&vsize);
	}

	time(&secondTime);

	//PetscPrintf(PETSC_COMM_WORLD,"End of Computation\n");

	PetscPrintf(PETSC_COMM_WORLD," TOTAL ELLAPSED TIME = %f seconds", difftime(secondTime,firstTime));


	/* end of the programm, deallocate structures and arrays */
 	MPI_Barrier(MPI_COMM_WORLD);
	VecDestroy(&v);
	MatDestroy(&A);

	/* free petsc arrays before finalizing slepc */
	PetscFree(com.in_received_buffer);
	PetscFree(com.out_sended_buffer);
	PetscFree(com.array_out_sended_buffer);
	PetscFree(com.array_in_received_buffer);

	/* free communicators and finalize */
	MPI_Barrier(MPI_COMM_WORLD); //wait for all
	SlepcFinalize(); //finalize petsc
	/* destroy communicators */
	mpi_lsa_com_free(&com);

	MPI_Finalize(); //finalize

	return 0;
}
