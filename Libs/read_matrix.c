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

 #include "read_matrix.h"

PetscErrorCode read_matrix_vector(Mat * A, Vec * v, int * communicator){
	char filea[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flaga,flagb;
	PetscViewer fd;
	PetscInt size,sizea;
	PetscScalar scal;

	ierr=PetscOptionsGetString(PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);
	if (!flaga) {
		sprintf(err,"Error : mfile is not properly set -> %s\n",filea);
		SETERRQ(*communicator,(PetscErrorCode) 83,err);
	}
	/* read matrix file */
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
///	PetscPrintf(PETSC_COMM_WORLD,"MatLoad gone OK\n");
	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",size,sizea);

	ierr=PetscOptionsGetString(PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);

	if (!flagb) {
		/* the user did not provide a vector, so generate it*/
		generate_random_seed_vector(size, -10.0, 10.0, 0, v);
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
