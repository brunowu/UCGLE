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
	
	mode_t process_mask = umask(0);
	/* init of MPI and MPI Utils */
	MPI_Init(&argc,&argv);
	mpi_lsa_init(argc,argv,&com);

	MPI_Barrier(MPI_COMM_WORLD);

	/* Launch PETSc */
	PETSC_COMM_WORLD=com.com_group; // give group communicator as main petsc communicator
	ierr=SlepcInitialize(&argc,&argv,(char *)0,help);	CHKERRQ(ierr);

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
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_nopc",&nols,&flag);CHKERRQ(ierr);
	if(!flag) nols=1;
	else if(nols==0) PetscPrintf(com.com_world,"]> not using LSA preconditionning %d\n",nols);

	if(com.color_group==GMRES_GR){
	                PetscPrintf(com.com_group,"]> Launched GMRES\n");
		launchGMRES(&com, &v, &A);
	//	PetscPrintf(com.com_group,"]> Launched GMRES\n");

	}else if(com.color_group==ARNOLDI_GR && nols!=0){
		Arnoldi(&com,&A,&v)	;
	}	else if(com.color_group==FATHER_GR && nols!=0){
		Father(&com,&v);
	} else if(com.color_group==LS_GR && nols!=0){
		VecGetSize(v,&vsize);
		LSQR(&com,&vsize);
	}

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
