static char help[] = "Solves a linear system in parallel with KSP.\n\n";

#include "slave.h"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );
  
  PetscErrorCode ierr;
  Vec x,b;
  Mat A;
  KSP ksp;
  int parentsize;
  int currentsize;
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD,"]> Initializing PETSc/SLEPc\n");

  int myrank;
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &currentsize );
  
  MPI_Comm parentcomm;
  MPI_Comm_get_parent( &parentcomm );
  MPI_Comm_remote_size(parentcomm, &parentsize);
  
  PetscPrintf(PETSC_COMM_WORLD,"VALDIATION ]> Parent size is %d\n", parentsize);
  MPI_Comm intra_comm;
  MPI_Intercomm_merge( parentcomm, 1, &intra_comm );

  ierr=loadInputs(&A,&b,&x);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"]> Data loaded\n");

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver settings done\n");

  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver Launching solving process\n");
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver System solved\n");

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"]> Cleaned structures, finalizing\n");

  int exit_type = 666;
  int i;
  MPI_Status status;
  int flg;
  int sended = 0;
  MPI_Request *type_request;
  type_request = malloc(parentsize*sizeof(MPI_Request));

  PetscPrintf(PETSC_COMM_WORLD,"ERROR CHECK >>> 1-SLAVE: my size is %d \n", currentsize);
  
  for(i = 0; i < parentsize;i++){
  	MPI_Isend(&exit_type,1, MPI_INT, i,i, intra_comm, &type_request[i]);
  	sended++;
  }

  PetscPrintf(PETSC_COMM_WORLD,"ERROR CHECK >>> 1-SLAVE: sended number to MASTER is %d \n",sended);
  
  PetscPrintf(PETSC_COMM_WORLD,"\nmsg>>> Child 1 send msg = %d \n", exit_type);
  //PetscFinalize();

  MPI_Comm_free(&parentcomm);
  MPI_Comm_free(&intra_comm);
  PetscFinalize();
  //MPI_Finalize();
//  printf("\nlllalalala\n\n");
  return 0;
}


