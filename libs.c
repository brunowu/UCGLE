#include "libs.h"

PetscErrorCode loadInputs(Mat * A, Vec * b, Vec * x){
  PetscErrorCode ierr;
  PetscInt sizex,sizey;
  char bfile[]="-bfile";
  char xfile[]="-xfile";
  
  //load data files
  ierr=loadMatrix(A);CHKERRQ(ierr);
  ierr=loadVector(bfile,b);CHKERRQ(ierr);
  if(*b==NULL) {
    PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector b\n");
    ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
    ierr=generateVectorRandom(sizex,b);CHKERRQ(ierr);
  }
  ierr=loadVector(xfile,x);CHKERRQ(ierr);
  if(*x==NULL) {
    PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector x\n");
    ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
    ierr=generateVectorRandom(sizex,x);CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode loadMatrix(Mat * A){
  char file[PETSC_MAX_PATH_LEN];
  char err[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  PetscBool flag;
  PetscViewer fd;
  PetscInt sizex,sizey;

  /*check args, if no matrix then no work... matrix file is mandatory*/
 
 ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD,"no file read \n");
    sprintf(err,"Error : mfile is not properly set -> %s\n",file);
    SETERRQ(PETSC_COMM_WORLD,(PetscErrorCode)83,err);
  }
  /* read matrix file */
  PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",file);

  ierr=MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
  ierr=MatSetType(*A,MATAIJ);CHKERRQ(ierr);

  ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr=MatLoad(*A,fd);CHKERRQ(ierr);
  ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",sizex,sizey);

  return 0;
}

PetscErrorCode loadVector(char * type_v,Vec * b){
  char file[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  PetscBool flag;
  PetscViewer fd;
  PetscInt size;

  // check if there is a vec file, vector is not mandatory
  ierr=PetscOptionsGetString(NULL,PETSC_NULL,type_v,file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD,"Error : %s is not properly set\n",type_v);
    *b = NULL;
  }else{
    PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",file);
    ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr=VecLoad(*b,fd);CHKERRQ(ierr);
    ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
    ierr=VecGetSize(*b,&size);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
  }

  return 0;
}

PetscErrorCode generateVectorRandom(int size, Vec * v){
  PetscErrorCode ierr;

  ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
  ierr=generateVector(size,v);CHKERRQ(ierr);
  ierr=VecSetRandom(*v,PETSC_NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Generated Random Vector of size : %d\n",size);

  return 0;
}

PetscErrorCode generateVectorNorm(int size, Vec * v){
  PetscScalar scal;
  PetscErrorCode ierr;

  ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
  ierr=generateVector(size,v);CHKERRQ(ierr);
  scal=1.0/size;
  ierr=VecSet(*v,scal);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Generated Norm Vector of size : %d\n",size);

  return 0;
}

PetscErrorCode generateVector(int size, Vec * v){
  PetscErrorCode ierr;

  ierr=VecCreate(PETSC_COMM_WORLD,v);CHKERRQ(ierr);
  ierr=VecSetSizes(*v,PETSC_DECIDE,size);CHKERRQ(ierr);
  ierr=VecSetFromOptions(*v);CHKERRQ(ierr);
  /*initiate the vector to its norm*/

  return 0;
}
