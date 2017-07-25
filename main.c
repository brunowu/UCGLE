/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "main.h"

static const char help[] = "Solve";

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char ** argv){
	com_lsa com;
	int vsize;
	Mat A,Af; 	
	Vec v,vf;
	PetscErrorCode ierr;
	PetscInt nols;
	PetscBool flag, gft_flg;
	int	non_lsa, size, rank;
	/* init of MPI and MPI Utils */
	MPI_Init(&argc,&argv);
  //      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	PetscPrintf("MPI WORLD SIZE = %d \n", size);
	non_lsa = mpi_lsa_init(argc,argv,&com);
	MPI_Barrier(MPI_COMM_WORLD);
	ierr=SlepcInitialize(&argc,&argv,(char *)0,help);       CHKERRQ(ierr);
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);	
        PetscPrintf(PETSC_COMM_WORLD,"MPI WORLD SIZE = %d \n", size);
	if(!non_lsa){

	/* Launch PETSc */
	PETSC_COMM_WORLD=com.com_group; // give group communicator as main petsc communicator
//	ierr=SlepcInitialize(&argc,&argv,(char *)0,help);	CHKERRQ(ierr);

//	PetscPrintf(com.com_world,"]> Initializing PETSc/SLEPc\n");

/*	VecCreate()*/

	/* programm options */
	ierr=PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the GMRES/LS-Arnoldi solver", "KSP");CHKERRQ(ierr);

	/* now read data files */
	MPI_Barrier(MPI_COMM_WORLD);
	PetscPrintf(com.com_world,"]> Opening input data files\n");
	MPI_Barrier(MPI_COMM_WORLD);
	
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

	mpi_lsa_com_vecsize_init(&com,&v);

	MPI_Barrier(MPI_COMM_WORLD);
	
	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_nopc",&nols,&flag);CHKERRQ(ierr);
	if(!flag) nols=1;
	else if(nols==0) PetscPrintf(com.com_world,"]> not using LSA preconditionning %d\n",nols);
	
	PetscOptionsHasName(NULL,NULL,"-GMRES_FT",&gft_flg);
	if(!gft_flg)
	{

		if(com.color_group==GMRES_GR){
			PetscPrintf(com.com_group,"]> Launching GMRES\n");
			launchGMRES(&com, &v, &A);
		}else if(com.color_group==ARNOLDI_GR && nols!=0){
			PetscPrintf(com.com_group,"]> Launching Arnoldi\n");
			Arnoldi(&com,&A,&v)	;
			}
		else if(com.color_group==FATHER_GR && nols!=0){
			PetscPrintf(com.com_group,"]> Launching Father\n");
			Father(&com,&v);
		} else if(com.color_group==LS_GR && nols!=0){
			VecGetSize(v,&vsize);
			PetscPrintf(com.com_group,"]> Launching LS\n");
			LSQR(&com,&vsize);
		}

	}
	else{


               if(com.color_group==GMRES_GR){
                        PetscPrintf(com.com_group,"]> Launching GMRES\n");
                 //       launchGMRES(&com, &v, &A);
                }else if(com.color_group==ARNOLDI_GR && nols!=0){
                        PetscPrintf(com.com_group,"]> Launching Arnoldi\n");
                        launchGMRES(&com, &v, &A);
                        }
                else if(com.color_group==FATHER_GR && nols!=0){
                        PetscPrintf(com.com_group,"]> Launching Father\n");
                        Father(&com,&v);
                } else if(com.color_group==LS_GR && nols!=0){
                        VecGetSize(v,&vsize);
                        PetscPrintf(com.com_group,"]> Launching LS\n");
                        LSQR(&com,&vsize);
                }


	}

        MPI_Barrier(MPI_COMM_WORLD);
        VecDestroy(&v);
        MatDestroy(&A);

        /* free petsc arrays before finalizing slepc */
        PetscFree(com.in_received_buffer);
        PetscFree(com.out_sended_buffer);
        PetscFree(com.array_out_sended_buffer);
        PetscFree(com.array_in_received_buffer);
        MPI_Barrier(MPI_COMM_WORLD); //wait for all
        SlepcFinalize(); //finalize petsc
        /* destroy communicators */
        mpi_lsa_com_free(&com);
	
	}
	else{
//        ierr=SlepcInitialize(&argc,&argv,(char *)0,help);CHKERRQ(ierr);
//        ierr=PetscInitialize(&argc,&argv,(char *)0,help);CHKERRQ(ierr);
	read_matrix_vector(&A, &v,PETSC_COMM_WORLD);
	classicalGMRES(&v, &A);
        MPI_Barrier(MPI_COMM_WORLD);
        VecDestroy(&v);
        MatDestroy(&A);
        MPI_Barrier(MPI_COMM_WORLD); //wait for all
        SlepcFinalize(); //finalize petsc

	}
	
	MPI_Finalize(); //finalize

	return 0;
}
