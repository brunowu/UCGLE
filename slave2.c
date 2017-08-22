static const char help[] = "Solve Non-Hermitian eigenvalues problem by Arnoldi, options array_in_received_buffer\n\
\t-mfile matrix_file (matrix file in PETSc bin format, this is mandatory)\n\
\t-xfile initial_guess_file (in PETSc bin format)\n";

#include <stdio.h> 
#include <mpi.h>
#include "slave2.h"

int main( int argc, char *argv[] ) {

  PetscErrorCode ierr;
  Vec x,b;
  Mat A;
  EPS eps;
  PetscInt its, nev, nconv;
  EPSType type;

  MPI_Init( &argc, &argv );
  ierr=SlepcInitialize(&argc,&argv,PETSC_NULL,help);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"\n\n]> Initializing SLEPc\n");
  int myrank;
  int currentsize;
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
//  MPI_Comm_size(MPI_COMM_WORLD, &currentsize);

  MPI_Comm parentcomm;
  MPI_Comm_get_parent( &parentcomm );

  MPI_Comm intra_comm2;
  MPI_Intercomm_merge( parentcomm, 1, &intra_comm2 );

  PetscPrintf(PETSC_COMM_WORLD, "VALIDATION ]> This is hello from SLAVE2\n");

  ierr=loadInputs(&A,&x,&b);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"]> Data loaded\n");

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, A, NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetType(eps,EPSARNOLDI);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver settings done\n");

  /*Solve the problem*/
  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver Launching solving process\n");

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver System solved\n");

  /*Get some informations of resolution*/

  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);
  /*Display the solution*/
  EPSGetConverged(eps,&nconv);
  PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);
  /*Clean*/
//  printf("\n\n HELLLLLLO \n\n");
  int exit_type = 0;
  MPI_Status status;
  int flag = 0,count;
  int end = 0;
  while(!end){
  	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,intra_comm2,&flag, &status);
  	if(flag){
  		MPI_Recv(&exit_type,1, MPI_INT, status.MPI_SOURCE,status.MPI_TAG,intra_comm2,&status);
        	printf("\nmsg>>> Child 2 recv msg = %d \n", exit_type);  
	}
	if(exit_type == 666){
		end = 1;
		break;
	}
  }
  MPI_Comm_size(MPI_COMM_WORLD, &currentsize);
  PetscPrintf(PETSC_COMM_WORLD,"ERROR CHECK >>> 2-SLAVE: my size is %d\n",currentsize);

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"]> SLEPC finalized\n");
 
  SlepcFinalize();
  MPI_Comm_free(&parentcomm);
  MPI_Comm_free(&intra_comm2);
 
    //printf("\n\n HELLLLLLO \n\n");
  //MPI_Finalize();
//  printf("\n\n HELLLLLLO \n\n");
  return 0;
}

