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

 #include "mpi_lsa_com.h"

/***********************************************************************************
Exchange information about the vector location within each intracommunicator to allow asynchronous communication.
************************************************************************************/
int mpi_lsa_com_vecsize_init(com_lsa * com, Vec * v){
	PetscErrorCode ierr;
	int * size1,*size2;
	int number1,number2,i;
	MPI_Comm com1,com2;
	PetscInt global_size;
	PetscInt local_size;
	int group_size;

	/* get vector size in global and local */
	ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);
	ierr=VecGetSize(*v,&global_size);CHKERRQ(ierr);

	/* get communicators sizes */
	MPI_Comm_remote_size(com->out_com,&(com->out_number));
	MPI_Comm_remote_size(com->in_com,&(com->in_number));

	MPI_Barrier(MPI_COMM_WORLD);
	group_size = com->size.com[com->color_group];
	/* alloc arrays depending on the number of remote processors */
	com->out_size=malloc(sizeof(int)*(com->out_number));
	com->in_size=malloc(sizeof(int)*(com->in_number));

	/* because we need to select the proper communicator between the two that link
	 * each intracom we will copy it into a variable */
	 //to make sure that com1 always size = 1

	if(com->color_group==GMRES_GR || com->color_group==ARNOLDI_GR){ // if many
		com1=com->in_com;
		com2=com->out_com;
		size1=com->in_size;
		size2=com->out_size;
		number1=com->in_number;
		number2=com->out_number;
	} else { // if one: ls or father processes
		com1=com->out_com;
		com2=com->in_com;
		size1=com->out_size;
		size2=com->in_size;
		number1=com->out_number;
		number2=com->in_number;
	}

	com->vec_requests=malloc(sizeof(MPI_Request)*(com->out_number));
	com->type_requests=malloc(sizeof(MPI_Request)*(com->out_number));
	com->array_requests=malloc(sizeof(MPI_Request)*(com->out_number));

	com->vec_in_disp=malloc(sizeof(int)*(com->in_number));
	com->vec_out_disp=malloc(sizeof(int)*(com->out_number));
	com->vec_color_disp=malloc(sizeof(int)*group_size);

	com->out_sended=0;
	com->in_received=0;

	PetscMalloc(local_size*sizeof(PetscScalar),&(com->out_sended_buffer));
	PetscMalloc(local_size*sizeof(PetscScalar),&(com->in_received_buffer));

	PetscMalloc(local_size*sizeof(PetscScalar),&(com->array_out_sended_buffer));
	PetscMalloc(local_size*sizeof(PetscScalar),&(com->array_in_received_buffer));

	com->out_vec_sended=0;
	com->array_out_sended=0;

	/* now process to the exchange */
	MPI_Allgather(&local_size,1,MPI_INT,size1,number1,MPI_INT,com1);
	MPI_Allgather(&local_size,1,MPI_INT,size2,number2,MPI_INT,com2);


	// displacements for each in or out node
	com->vec_in_disp[0]=0;
	for(i=1;i<com->in_number;i++){
		com->vec_in_disp[i]=com->vec_in_disp[i-1]+com->in_size[i-1];
	}
	com->vec_out_disp[0]=0;
	for(i=1;i<com->out_number;i++){
		com->vec_out_disp[i]=com->vec_out_disp[i-1]+com->out_size[i-1];
	}
	com->vec_color_disp[0]=0;

	return 0;
}

/***********************************************************************************
*	send data to one or many nodes through an intercommunicator
************************************************************************************/
int mpi_lsa_com_vec_send(com_lsa * com, Vec * v){
	int i,incr,tmp_int,vsize,send_size;
	PetscScalar * array;
	PetscScalar tmp_global,tmp_local;
	PetscErrorCode ierr;
	MPI_Status status;
	int flag;
	int out_com_size;
	Vec x;

    tmp_int=0;
	ierr=VecDuplicate(*v,&x);CHKERRQ(ierr);
	ierr=VecCopy(*v,x);CHKERRQ(ierr);

	/* check if previous requests where completed */
	for(i=0;i<com->out_vec_sended;i++){

	MPI_Test(&(com->vec_requests[i]),&flag,&status);
		/* if not cancel it */
		if(!flag){
			/* if no vector was sent by any of the group nodes*/
			tmp_local=(PetscScalar)com->out_vec_sended;
			tmp_global=(PetscScalar)0;

			/*now proceed to the communication*/
			/* we comput the sum of tmp local of the group within each node of the group*/
			MPI_Allreduce(&tmp_local,&tmp_global,1,MPIU_SCALAR,MPI_SUM,com->com_group);
			/* determines the size of com_group*/
			MPI_Comm_size(com->com_group,&vsize);

			/* we update the vector to send to the latest version */
			if(((int)(PetscRealPart(tmp_global)))==vsize){
				for(i=0;i<com->out_vec_sended;i++)
					MPI_Cancel(&(com->vec_requests[i]));
					com->out_vec_sended=0;
				}else{
					return 1;
			  }
			}else{
				com->out_vec_sended--;
			}
		}

	MPI_Barrier(PETSC_COMM_WORLD);

	/* extract the vector array */
	/* first we get the size of the array (stored in the local memory)
	 logiquement ce doit être un VecGetLocalSize */
	ierr=VecGetLocalSize(x,&vsize);CHKERRQ(ierr);
	ierr=PetscMalloc(sizeof(PetscScalar) * vsize, &array);CHKERRQ(ierr);
	ierr=VecGetArray(x,&array);CHKERRQ(ierr);

	/* we copy the array extracted to the sended_buffer */
	for(i=0;i<vsize;i++)
		(com->array_out_sended_buffer)[i]=(PetscScalar)array[i];

	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		/* compute array displacement */
		if(i<1)
			incr=0;
		else
			incr=com->vec_out_disp[i];
			MPI_Comm_size(com->com_group,&tmp_int);

		send_size=com->out_size[i];
		if(vsize<com->out_size[i]) send_size=vsize;

		/* we send a portion of data */
		MPI_Isend(&(com->array_out_sended_buffer)[incr],send_size,MPIU_SCALAR,i,com->rank_group,com->out_com,&(com->vec_requests[i]));
		com->out_vec_sended++;
	}

	/* restore array */
	ierr=VecRestoreArray(x,&array);CHKERRQ(ierr);
	return 0;
}

/* for the validatig funcion check on the bottom */
int mpi_lsa_com_vec_recv(com_lsa * com, Vec * v){
	int flag,count,i;
	PetscInt size,tmp_int;
	PetscErrorCode ierr;
	MPI_Status status,statas;

	MPI_Comm_size(com->com_group,&tmp_int);
	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);
	/* did we received something ? */
	if(!flag){
		return 1; // didn't received anything
	}
	/* just verify coerancy of data */
	MPI_Get_count(&status,MPIU_SCALAR,&count);
	if(count == 1 || count ==0)
		return 1;
	ierr=VecGetSize(*v,&size);CHKERRQ(ierr);

	/* get the array, must be done because we can't restore the array otherwise
	 * this operation is faster than setting the values of the array (i think).
	 * the array pointer is only pointing to the local vector memory of the node */
	MPI_Recv(&(com->array_in_received_buffer[com->vec_in_disp[status.MPI_SOURCE]]),count,MPIU_SCALAR,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&statas);
	com->in_received++;

	for( i = 0; i < count; ++i){
		ierr=VecSetValue(*v, (PetscInt)i+com->vec_in_disp[status.MPI_SOURCE],com->array_in_received_buffer[i], INSERT_VALUES);CHKERRQ(ierr);
/*		printf("com->in_received_buffer[%d] = %f\n", i,com->array_in_received_buffer[i]);*/
	}
	return 0;
}

int mpi_lsa_com_array_send(com_lsa * com, int * size, PetscScalar * data){
	MPI_Status status;
	int flag,i,vsize;
	PetscScalar tmp_global,tmp_local;

	/* check if previous requests where completed */
	for(i=0;i<com->array_out_sended;i++){
		MPI_Test(&com->array_requests[i],&flag,&status);
		/* if not cancel it */
		if(!flag){
			/* if no vector was sent by any of the group nodes*/
			tmp_local=(PetscScalar)com->array_out_sended;
			tmp_global=(PetscScalar)0;
			MPI_Allreduce(&tmp_local,&tmp_global,1,MPIU_SCALAR,MPI_SUM,com->com_group);

			MPI_Comm_size(com->com_group,&vsize);

			/* we update the vector to send to the latest version */
			if(((int)(PetscRealPart(tmp_global)))==vsize){
				for(i=0;i<com->array_out_sended;i++){
					MPI_Cancel((MPI_Request*)&(com->array_requests[i]));
				}
				com->array_out_sended=0;
			}else{
				return 1;
			}
		}else{
			com->array_out_sended--;
		}
	}
	/*fill the send buffer with the data */
	for(i=0;i<*size;i++){
		com->array_out_sended_buffer[i]=(PetscScalar)data[i];
	}
	com->nbr_array_sended +=1;
	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		/* we send a portion of data */
	       
 		printf("Send size=%d (%d scalars) to %d \n",*size,*size/8,i);

		MPI_Isend(com->array_out_sended_buffer,*size,MPIU_SCALAR,i,com->nbr_array_sended + i,com->out_com,&com->array_requests[i]);
		com->out_vec_sended++;
	}
	return 0;
}

int mpi_lsa_com_array_recv(com_lsa * com, int * size, PetscScalar * data){
	MPI_Status status;
	int flag;

	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);
	MPI_Get_count(&status,MPI_INT,size);
	/* did we received something ? */
	if(!flag || *size==1)
		return 1; // didn't received anything

	/* how large will be the array */
	MPI_Get_count(&status,MPIU_SCALAR,size);

	printf("%d Receive size=%d (%d scalars) msg Number = %d\n", com->rank_world,*size,*size/8, status.MPI_TAG-com->rank_group);

	/* receive the data array */
	MPI_Recv(data,*size,MPIU_SCALAR,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&status);
	return 0;
}

int mpi_lsa_send_vec(com_lsa * com, Vec * x){
  PetscScalar * tmp_array;
  int i, master_rank, my_rank, group_id, local_size, global_size, ierr, group_size, flag;
  MPI_Status status;
  global_size = 0;
  Vec * v;

  ierr = VecDuplicate(*x, v);CHKERRQ(ierr);
  MPI_Comm_rank(com->com_group, &my_rank);
  MPI_Comm_size(com->com_group, &group_size);
  group_id = com->color_group;
  master_rank = com->master.com[group_id];

  ierr = VecGetLocalSize(*v, &local_size);CHKERRQ(ierr);
  tmp_array = malloc(local_size * sizeof(PetscScalar));
  ierr = VecGetArray(*v, &tmp_array); CHKERRQ(ierr); // good

  if(com->rank_world == master_rank){
		 i = 0;
		 while(com->out_vec_sended > 0){
				i = com->out_vec_sended;
				MPI_Test(&(com->vec_requests[i]),&flag,&status);
				/* if not cancel it */
				if(!flag){
					/* we update the vector to send to the latest version */
					MPI_Cancel(&(com->vec_requests[i]));
					com->out_vec_sended--;
				}
		}
		global_size = 0;
    for(i = 0; i < group_size; ++i){
        global_size += com->vec_color_sizes[i];
		}
    com->out_sended_buffer = malloc(global_size * sizeof(PetscScalar));
  }

  /* first we gather all the vecs into one in master*/
  if(group_size > 1){
    MPI_Gather(tmp_array, local_size, MPIU_SCALAR, com->out_sended_buffer, global_size,MPIU_SCALAR, 0, com->com_group);
	}else
		com->out_sended_buffer = tmp_array;

  /* Now the root(master) will proceed to the sending*/
  if(com->rank_world == master_rank) {
    if(com->color_group < 3)
      MPI_Isend(com->out_sended_buffer, global_size, MPIU_SCALAR, com->master.com[group_id + 1], com->rank_world,  MPI_COMM_WORLD, &(com->vec_requests[0]));
    else
      MPI_Isend(com->out_sended_buffer, global_size, MPIU_SCALAR, com->master.com[0], com->rank_world, MPI_COMM_WORLD, &(com->vec_requests[0]));
  }
	ierr = VecRestoreArray(*v, &tmp_array);CHKERRQ(ierr);
	return 0;
}

int mpi_lsa_receive_vec(com_lsa * com, Vec * v){

	PetscScalar * tmp_array;
	int i, master_rank, my_rank, group_id, local_size, global_size, ierr, group_size, flag;
  MPI_Status status, statas;
  int * indices;

  global_size = 0;
  MPI_Comm_rank(com->com_group, &my_rank);
  MPI_Comm_size(com->com_group, &group_size);
  group_id = com->color_group;
  master_rank = com->master.com[group_id];

	local_size = com->vec_color_sizes[my_rank];
  tmp_array = malloc(local_size * sizeof(PetscScalar));
	if(group_size > 1){
	 	indices = malloc(group_size * sizeof(int));
		indices[0] = 0;
		for(i = 1; i < group_size; ++i){
		  indices[i] = indices[i-1] + (uintptr_t)com->vec_color_sizes;
		}
	}
	com->recv_flag = -1;
	if( (com->rank_world == com->master.com[group_id]) && (my_rank == 0) ){

	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
	MPI_Get_count(&status,MPIU_INT,&global_size);

	/* did we received something ? */
	if(!flag || global_size==1){
		com->recv_flag = 0;
		MPI_Bcast(&(com->recv_flag), 1, MPI_INT, 0, com->com_group);
		return 1; // didn't received anything
	}
	com->recv_flag = 1;
	MPI_Bcast(&(com->recv_flag), 1, MPI_INT, 0, com->com_group);

	/* we do it for one side because of LS that receives different sizes and also to be sure that the size is correct*/
	com->in_received_buffer = malloc(global_size * sizeof(PetscScalar));
	MPI_Recv(com->in_received_buffer,global_size,MPIU_SCALAR,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD, &statas);
	com->in_received++;

	}else{
	 	MPI_Bcast(&(com->recv_flag), 1, MPI_INT, 0, com->com_group);
	 	if(!flag ){
			return 1; // didn't received anything
		}
	}
	MPI_Barrier(com->com_group);
	if(group_size > 1){
		MPI_Scatterv( com->in_received_buffer, com->vec_color_sizes, indices, MPIU_SCALAR,
		tmp_array, local_size, MPIU_SCALAR, 0, com->com_group);
	}
	return 0;
}

int mpi_lsa_com_free(com_lsa * com){
	/* free arrays */
	free(com->out_size);
	free(com->in_size);
	free(com->vec_requests);
	free(com->type_requests);
	free(com->array_requests);
	free(com->vec_in_disp);
	free(com->vec_out_disp);
	/* and mpi communicators */
	MPI_Comm_free(&(com->out_com));
	MPI_Comm_free(&(com->in_com));
	MPI_Comm_free(&(com->com_group));

	return 0;
}

int mpi_lsa_com_type_recv(com_lsa * com, int * type){
	MPI_Status status;
	int flag,size;

	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);

	if(!flag )
		return 1; // didn't received anything
		MPI_Get_count(&status,MPI_INT,&size);

	if(size!=1)
		return 1;
		MPI_Recv(type,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&status);

	return 0;
}

int mpi_lsa_com_type_send(com_lsa * com, int * type){
	MPI_Status status;
	int flag,i;

	/*check if previous requests where completed */
	for(i=0;i<com->out_sended;i++){
		MPI_Test(&(com->type_requests[i]),&flag,&status);
		/* if not cancel it */
		if(!flag)
			MPI_Cancel(&(com->type_requests[i]));
	}
	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		MPI_Isend(type,1,MPI_INT,i,i,com->out_com,&(com->type_requests[i]));
		com->out_sended++;
	}
	return 0;
}


int mpi_lsa_com_vec_recv_validate(com_lsa * com, Vec * v, int size){
	int vector_chunks;
	int my_com_size;
	PetscScalar  *array;
	int global_size,local_size,i;
	PetscErrorCode ierr;

	/* first we get the number of vector portions that the group should receive */
	MPI_Comm_size(com->com_group,&my_com_size);

	if(my_com_size == 1){ //this means that the process is alone in his group
	  ierr=VecGetSize(*v,&local_size);CHKERRQ(ierr);
	  if(size == local_size){
	   	ierr=VecGetArray(*v,&array);CHKERRQ(ierr);
		for(i=0;i<local_size;i++)
			array[i]=com->in_received_buffer[i];
		ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
		com->in_received=0;
		return 0;
	  } else return 1;
	}

	if(my_com_size > 1){
	  global_size=0;

	  ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);
	  MPI_Allreduce( &local_size, &global_size, 1, MPI_INT, MPI_SUM, com->com_group);
	  if(global_size == size){
	   	ierr=VecGetSize(*v,&local_size);CHKERRQ(ierr);
		  if(size == local_size){
		    ierr=VecGetArray(*v,&array);CHKERRQ(ierr);
			  for(i=0;i<local_size;i++)
					array[i]=com->in_received_buffer[i];
			  ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
			  com->in_received=0;
		   return 0;
		  }
		}else return 1;
	}
	return 1; //this happends when my_com_size has strenge value "less than 0"
}

int setting_out_vec_sizes( com_lsa * com, Vec * v){
	int i, tag = 1707, flag = 0, ierr;
	PetscInt local_size;
	int master_rank, group_id, my_rank, group_rank;
	MPI_Request * request;
	MPI_Status  status;

	request = malloc( 2 * sizeof(MPI_Request));
	group_id = com->color_group;
	master_rank = com->master.com[group_id];
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_rank(com->com_group,&group_rank);
	com->vec_color_sizes = malloc( com->size.com[group_id] * sizeof(int));
	com->vec_out_sizes = malloc(com->out_number * sizeof(int));

	ierr=VecGetSize(*v, &(com->vec_color_sizes[group_rank]));CHKERRQ(ierr);
	if(com->size.com[group_id] > 1){
		MPI_Allgather(&(com->vec_color_sizes[group_rank]), 1, MPIU_INT,
                  com->vec_color_sizes, com->size.com[group_id], MPIU_INT,
                  com->com_group);
	}
  if(my_rank == master_rank){
   	if(com->color_group != 0)
   		MPI_Isend(com->vec_color_sizes, com->size.com[group_id], MPI_INT,
   		          com->master.com[group_id - 1], tag + my_rank, MPI_COMM_WORLD, &request[0]);
    else
   		MPI_Isend(com->vec_color_sizes, com->size.com[group_id], MPI_INT,
   		          com->master.com[3], tag + my_rank, MPI_COMM_WORLD, &request[0]);
   	if(com->color_group != 3)
   		MPI_Recv(com->vec_out_sizes, com->out_number, MPI_INT,
   		         com->master.com[group_id + 1], tag + com->master.com[group_id + 1], MPI_COMM_WORLD, &status);
   	else
   		MPI_Recv(com->vec_out_sizes, com->out_number, MPI_INT,
   		         com->master.com[0], tag , MPI_COMM_WORLD, &status);
  }
	MPI_Ibcast(com->vec_out_sizes, com->out_number, MPI_INT, 0, com->com_group, &request[1]);
	return 0;
}
