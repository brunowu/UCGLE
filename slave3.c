#include "slave3.h"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );

  int myrank;
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Comm parentcomm;
  MPI_Comm_get_parent( &parentcomm );

  MPI_Comm intra_comm3;
  MPI_Intercomm_merge( parentcomm, 1, &intra_comm3 );

  if ( myrank == 0 ) {
    int tmp_size;
    MPI_Comm_remote_size( parentcomm, &tmp_size );
    //  printf( "child size of parent comm world = %d\n", tmp_size );

    MPI_Comm_size( MPI_COMM_WORLD, &tmp_size );
    //  printf( "child size of child comm world = %d\n", tmp_size );

    MPI_Comm_size( intra_comm3, &tmp_size );
    //   printf( "child size of intra comm world = %d\n", tmp_size );
  }
  printf("VALIDATION ]> This is HELLO from SLAVE 3!!!\n");
  //  printf("I am rank %d in intra_rank \n", intra_rank);

  int exit_type = 0;
  MPI_Status status;
  int flag = 0,count;
  int end = 0;

  while(!end){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,intra_comm3,&flag, &status);
    if(flag){
      MPI_Recv(&exit_type,1, MPI_INT, status.MPI_SOURCE,status.MPI_TAG,intra_comm3,&status);
      printf("msg>>> Child 3 recv msg = %d \n", exit_type);
    }
    if(exit_type == 666){
      end = 1;
      break;
    }
  }

  MPI_Comm_free(&parentcomm);
  MPI_Comm_free(&intra_comm3);
  //MPI_Finalize();

  return 0;

}
