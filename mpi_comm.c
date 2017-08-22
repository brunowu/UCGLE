#include "mpi_comm.h"

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
