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


#include "arnoldi.h"

/* Compute cyclicly eigenvalue */
PetscErrorCode Arnoldi(com_lsa * com, Mat * A, Vec  *v){

	EPS eps; /* eigensolver context */
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	PetscInt end,first,validated;
	PetscErrorCode ierr;
	/* eigenvalues number is set to 100, can be changed if needed
	   we choosed to fix it because mallocs weren't working properly */
	PetscScalar eigenvalues[500], ei, er;
	PetscReal re,im,vnorm;
	PetscInt eigen_nb,j,i,size, taille;
	uintptr_t one=1;
	Vec initialv,nullv,*vs;
	PetscBool flag,data_load,data_export,continuous_export,load_any;
	int exit_type=0, counter = 0, l;
	int sos_type = 911;
	Vec vecteur_initial;
	PetscViewer viewer;

	char  ls_load_path[PETSC_MAX_PATH_LEN];
	PetscBool ls_load, ls_load_any;

	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_load",ls_load_path,PETSC_MAX_PATH_LEN,&ls_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_load_any",&ls_load_any);CHKERRQ(ierr);

	if(ls_load&&ls_load_any){
	  ls_load_any=PETSC_FALSE;
	  ls_load=PETSC_TRUE;
	}
	PetscBool need_new_init = PETSC_FALSE, exit = PETSC_FALSE;

	end=0;
	first=1;
	validated=1;

 if(!(ls_load^=ls_load_any)){

	sprintf(load_path,"./arnoldi.bin");
	sprintf(export_path,"./arnoldi.bin");

	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);

	/* create the eigensolver */
	ierr=EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	/* set the matrix operator */
	ierr=EPSSetOperators(eps,*A,PETSC_NULL);
	/* set options */
	ierr=EPSSetType(eps,EPSARNOLDI);
	ierr=EPSSetFromOptions(eps);CHKERRQ(ierr);
	/* duplicate vector properties */
	ierr=VecDuplicate(*v,&initialv);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&nullv);CHKERRQ(ierr);
	ierr=VecSet(nullv,(PetscScalar)0.0);CHKERRQ(ierr);

	ierr=VecSetRandom(initialv,PETSC_NULL);//initialize initial vector to random
	ierr=VecGetSize(initialv,&size);CHKERRQ(ierr);

	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen",&eigen_nb,&flag);CHKERRQ(ierr);
	if(!flag) eigen_nb=EIGEN_ALL;
	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_arnoldi_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_arnoldi_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);

	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_arnoldi_load_any",&load_any);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_arnoldi_cexport",&continuous_export);CHKERRQ(ierr);

	vs=malloc(size*sizeof(Vec));
	for(i=0;i<size;i++){
		ierr=VecDuplicate(*v,&vs[i]);CHKERRQ(ierr);
	}
	ierr=VecDuplicate(initialv,&vecteur_initial);CHKERRQ(ierr);

	while(!end){
		/*check if the program need to exit */
		if(exit == PETSC_TRUE)
			break;
		/* check if we received an exit message from Father*/
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
		/* check if we received some data from GMRES		*/
		if(!mpi_lsa_com_vec_recv(com, &initialv)){
				VecGetSize(initialv, &taille);
		}

		for(j=0;j<eigen_nb;j++){
			eigenvalues[j]=(PetscScalar)0.0;
		}

		if(data_load&&load_any){
		  load_any=PETSC_FALSE;
		  data_load=PETSC_TRUE;
		}
		ierr = VecAssemblyBegin(initialv);CHKERRQ(ierr);
  	ierr = VecAssemblyEnd(initialv);CHKERRQ(ierr);

		if(!(data_load^=load_any)){
		  ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);
		} else {
		  ierr=readBinaryVecArray(load_path,(int*)one,&initialv);CHKERRQ(ierr);
			data_load=PETSC_FALSE;
			load_any=PETSC_FALSE;
			ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);
		}
		ierr=EPSSolve(eps);CHKERRQ(ierr);
		/*construct new initial vector*/
		ierr=EPSGetInvariantSubspace(eps, vs);CHKERRQ(ierr);
		++counter;
		/* get the number of guessed eigenvalues */
		ierr=EPSGetConverged(eps,&eigen_nb);CHKERRQ(ierr);
		/* send them */
		for(j=0;j<eigen_nb;j++){
			ierr = EPSGetEigenvalue(eps,j,&er,&ei);CHKERRQ(ierr);
			#ifdef PETSC_USE_COMPLEX
			  re=PetscRealPart(er);
			  im=PetscImaginaryPart(er);
			#else
			  re=er;
			  im=ei;
			#endif
			eigenvalues[j]=(PetscScalar)re+PETSC_i*(PetscScalar)im;
		}
		if( eigen_nb != 0){
			mpi_lsa_com_array_send(com, &eigen_nb, eigenvalues);
			/*construct new initial vector*/
			ierr=VecCopy(vs[0],initialv);CHKERRQ(ierr);

			for(j=1;j<eigen_nb;j++){
				ierr=VecAYPX(initialv,(PetscScalar)1.0,vs[j]);
			}

			ierr=VecNorm(initialv,NORM_2,&vnorm);CHKERRQ(ierr);
			ierr=VecAYPX(initialv,(PetscScalar)(1.0/vnorm),nullv);CHKERRQ(ierr);

	  	if(continuous_export){
		  	ierr=writeBinaryVecArray(data_export?export_path:"./arnoldi.bin", 1, &initialv);
			}
			if(!mpi_lsa_com_type_recv(com,&exit_type)){
   			if(exit_type==666){
		 			end=1;
		 			mpi_lsa_com_type_send(com,&exit_type);
		 			need_new_init = PETSC_FALSE;
		 			exit = PETSC_TRUE;
					break;
  				}
			}
  	}else{
  		need_new_init = PETSC_TRUE;
 			while(need_new_init){
				 if(!mpi_lsa_com_vec_recv(com, &initialv)){
					 need_new_init = PETSC_FALSE;
		     }else{
		       if(!mpi_lsa_com_type_recv(com,&exit_type)){
			     		if(exit_type==666){
				    		end=1;
				    		mpi_lsa_com_type_send(com,&exit_type);
				    		need_new_init = PETSC_FALSE;
				    		exit = PETSC_TRUE;
				   			break;
				  		}
					 }
      	}
		  }
		  if(exit == PETSC_TRUE)
		     break;
  	}
	}

if(data_export){
	ierr=writeBinaryVecArray(export_path, 1, &initialv);
}
for(i=0;i<eigen_nb;i++){
	ierr=VecDestroy(&(vs[i]));CHKERRQ(ierr);
}
/* and destroy the eps */
ierr=EPSDestroy(&eps);CHKERRQ(ierr);
ierr=VecDestroy(&initialv);CHKERRQ(ierr);
ierr=VecDestroy(&nullv);CHKERRQ(ierr);
}

else{
	while(!end){
		/*check if the program need to exit */
		if(exit == PETSC_TRUE)
			break;
		/* check if we received an exit message from Father*/
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
   }
	}


return 0;
}
