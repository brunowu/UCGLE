static char help[] = "Solves a linear system in parallel with KSP.\n\n";

#include "slave.h"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );
  
  PetscErrorCode ierr;
  Vec x,b;
  Mat A;
  KSP ksp;

  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD,"]> Initializing PETSc/SLEPc\n");

  int myrank;
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Comm parentcomm;
  MPI_Comm_get_parent( &parentcomm );

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

  MPI_Comm_free(&parentcomm);
  MPI_Comm_free(&intra_comm);

  PetscFinalize();
  MPI_Finalize();

  return 0;
}


