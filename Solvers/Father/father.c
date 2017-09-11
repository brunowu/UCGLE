/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
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

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_ls_load",ls_load_path,PETSC_MAX_PATH_LEN,&ls_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_ls_load_any",&ls_load_any);CHKERRQ(ierr);

	if(ls_load&&ls_load_any){
	  ls_load_any=PETSC_FALSE;
	  ls_load=PETSC_TRUE;
	}

	ierr=VecDuplicate(*v,&vec_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&vec_tmp_receive);CHKERRQ(ierr);

	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
//		PetscPrintf(PETSC_COMM_WORLD,"\nFather exit type = %d \n", exit_type);
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

                  else if(exit_type==999){
	          	PetscPrintf(PETSC_COMM_WORLD, "\nFather receive and send the GMRES FT simulation msg\n");
                        end = 1;
			for(i=0;i<com->out_vec_sended;i++){

                                        MPI_Test(&(com->vec_requests[i]),&flag,&status);
                             
                                        if(!flag){
                                                MPI_Cancel(&(com->vec_requests[i]));
                                        }
                    }
                    mpi_lsa_com_type_send(com,&exit_type);
                        break;
                  }
		  else if(exit_type != 0){
		       //PetscPrintf(PETSC_COMM_WORLD, "\nGMRES RESOLUTION COUNT = %d\n",exit_type);
                        for(i=0;i<com->out_vec_sended;i++){

                                        MPI_Test(&(com->vec_requests[i]),&flag,&status);

                                        if(!flag){
                                                MPI_Cancel(&(com->vec_requests[i]));
                                        }
                    }

		       mpi_lsa_com_type_send(com,&exit_type);
//		       PetscPrintf(PETSC_COMM_WORLD, "\nGMRES has sent RESOLUTION COUNT = %d\n",exit_type);
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
