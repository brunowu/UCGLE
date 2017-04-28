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

#include "father.h"

/* Compute cyclicly eigenvalues, we create and destroy the EPS context inside the loop.
   This is not the best thing to do because it may decrease performances due to setup overhead */
PetscErrorCode Father(com_lsa * com, Vec * v){
	Vec vec_tmp_receive, vec_tmp;
	PetscInt end;
	PetscErrorCode ierr;
	int exit_type=0, i ;
	end=0;
	MPI_Status status;
	int flag;
	Vec vecteur_initial;
	int taille = 1;

	char  ls_load_path[PETSC_MAX_PATH_LEN];
	PetscBool ls_load, ls_load_any;

	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_load",ls_load_path,PETSC_MAX_PATH_LEN,&ls_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_load_any",&ls_load_any);CHKERRQ(ierr);

	if(ls_load&&ls_load_any){
	  ls_load_any=PETSC_FALSE;
	  ls_load=PETSC_TRUE;
	}

	ierr=VecDuplicate(*v,&vec_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&vec_tmp_receive);CHKERRQ(ierr);

	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;

		    for(i=0;i<com->out_vec_sended;i++){

					MPI_Test(&(com->vec_requests[i]),&flag,&status);
					/* if not cancel it */
					if(!flag){
					/* we update the vector to send to the latest version */
				 		MPI_Cancel(&(com->vec_requests[i]));
					}
		    }
		    mpi_lsa_com_type_send(com,&exit_type);
		  	break;
		  }
		}
		/* check if there's an incoming message */
		if(!mpi_lsa_com_vec_recv(com, &vec_tmp_receive) && !(ls_load^=ls_load_any)){
			if(com->in_received == com->in_number){
				ierr=VecCopy(vec_tmp_receive,vec_tmp);CHKERRQ(ierr);
				mpi_lsa_com_vec_send(com,&vec_tmp);
				com->in_received = 0;
			}
		}
	}
	return 0;
}
