/*
 This file is part of software for the implementation of UCGLE method, under the supervision of Serge G. Petiton
 <serge.petiton@univ-lille1.fr>.
 
 Copyright (C) 2011—. Pierre-Yves AQUILANTI and Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr> in Maison de la Simulation. 
 All rights reserved.
 
 Permission to use, copy, modify and distribute this software for personal and educational use is hereby granted
 without fee, provided that the above copyright notice appears in all copies and that both that copyright notice 
 and this permission notice appear in supporting documentation, and that the names of all authors are not used in 
 advertising or publicity pertaining to distribution of the software without specific, written prior permission. 
 Xinzhe WU and the author make no representations about the suitability of this software for any purpose. It is 
 provided "as is" without express or implied warranty.
 
 You should have received a copy of the GNU Lesser General Public License along with UCGLE.  If not, see 
 <http://www.gnu.org/licenses/>.

 For more information, contact with Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr>.
 
 */

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
