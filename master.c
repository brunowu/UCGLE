#include <stdio.h>
#include <mpi.h>

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  MPI_Comm child_comm, child_comm2;
  int  num_processes_to_spawn = 2;
  MPI_Comm_spawn( "./worker", argv+1,
		  num_processes_to_spawn, MPI_INFO_NULL,
		  0, MPI_COMM_WORLD,
		  &child_comm, MPI_ERRCODES_IGNORE );

  MPI_Comm_spawn( "./worker2", argv+1,
                  num_processes_to_spawn, MPI_INFO_NULL,
                  0, MPI_COMM_WORLD,
                  &child_comm2, MPI_ERRCODES_IGNORE );

  MPI_Comm intra_comm, intra_comm2;
  MPI_Intercomm_merge( child_comm, 0, &intra_comm );
  MPI_Intercomm_merge( child_comm2, 0, &intra_comm2 );

  int intra_rank, intra_rank2;
  MPI_Comm_rank(intra_comm, &intra_rank);
  MPI_Comm_rank(intra_comm2, &intra_rank2);

  printf("intra rank = %d \n\n\n",intra_rank);
  if ( rank == 0 ) {
    int tmp_size;
    int tmp_size2;
    MPI_Comm_size( intra_comm, &tmp_size );
    MPI_Comm_size( intra_comm2, &tmp_size2 );

    MPI_Comm_remote_size( child_comm, &tmp_size );
    MPI_Comm_remote_size( child_comm2, &tmp_size2 );

    MPI_Comm_size( MPI_COMM_WORLD, &tmp_size );
    MPI_Comm_size( MPI_COMM_WORLD, &tmp_size2 );

  }

  int exit_type = 0;
  MPI_Status status;
  int flag = 0, count;
  int end = 0;
  MPI_Request request; 
  while(!end){
  	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,intra_comm,&flag, &status);
  	if(flag){
  		MPI_Recv(&exit_type,1, MPI_INT, status.MPI_SOURCE,status.MPI_TAG,intra_comm,&status);
        	printf("\n### Master recv msg = %d \n", exit_type);  
	}
  	if(exit_type == 666){
		end = 1;
		MPI_Isend(&exit_type,1,MPI_INT, 1,0,intra_comm2,&request);
		break;
  	}
  }  
  MPI_Comm_free(&child_comm);
  MPI_Comm_free(&child_comm2);
  MPI_Comm_free(&intra_comm);
  MPI_Comm_free(&intra_comm2);
//  printf("\n\nhey\n\n");
  //MPI_Finalize();
  printf("\n\nhey\n\n");
  return 0;
}
