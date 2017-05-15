/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "data_rw.h"

PetscErrorCode writeBinaryScalarArray(const char * name, int nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;

  ierr=PetscBinaryOpen(name,FILE_MODE_WRITE,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedWrite(PETSC_COMM_WORLD,file_descriptor,array,nb,PETSC_SCALAR,PETSC_FALSE);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode writeBinaryVecArray(const char * name, int nb, Vec * array){
  PetscInt size,i;
  PetscScalar * tmp_array;
  PetscErrorCode ierr;

  ierr=VecGetSize(array[0],&size);CHKERRQ(ierr);

  for(i=0;i<nb;i++)
    ierr=VecGetArray(array[i],&tmp_array);CHKERRQ(ierr);
  ierr=writeBinaryScalarArray(name,nb*size,tmp_array);CHKERRQ(ierr);


  for(i=0;i<nb;i++)
    ierr=VecRestoreArray(array[i],&tmp_array);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode readBinaryScalarArray(const char * name, int * nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;
  long size;

  getFileSize(name,&size);

  if(*nb<=0) *nb=(int)size/((int)sizeof(PetscScalar));
  if(size/sizeof(PetscScalar)!=*nb) {
    return 1;
  }


  ierr=PetscBinaryOpen(name,FILE_MODE_READ,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedRead(PETSC_COMM_WORLD,file_descriptor,array,*nb,PETSC_SCALAR);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode readBinaryVecArray(const char * name, int * nb, Vec * array){
  PetscInt size,i;
  PetscScalar * tmp_array;
  PetscErrorCode ierr;
  ierr=VecGetSize(array[0],&size);CHKERRQ(ierr);

  for(i=0;i<*nb;i++)
    ierr=VecGetArray(array[i],&tmp_array);CHKERRQ(ierr);

  size*=(*nb);
  ierr=readBinaryScalarArray(name,&size,tmp_array);CHKERRQ(ierr);

  for(i=0;i<*nb;i++)
    ierr=VecRestoreArray(array[i],&tmp_array);CHKERRQ(ierr);

  return ierr;
}


PetscErrorCode getFileSize(const char * name, long * size){
  FILE * fptr;
  *size = 0L;

#ifdef LINUX
  struct stat fs;

  if(stat(name,&fs)!=0){
    perror("Cannot state file\n");
  }
  *size=fs.st_size;

#else
fptr=fopen(name,"rb");
  if(fptr!=NULL){
    fseek(fptr,0L,SEEK_END);
    *size = ftell(fptr);
    fclose(fptr);
  }
#endif

  return 0;
}
