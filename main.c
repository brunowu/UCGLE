/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */

#include "main.h"

static const char help[] = "Solve";

#undef __FUNCT__
#define __FUNCT__ "main"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  MPI_Comm child_comm, child_comm2;
  int  num_processes_to_spawn = 2;
  MPI_Comm_spawn( "./GMRES", argv+1,
		  num_processes_to_spawn, MPI_INFO_NULL,
		  0, MPI_COMM_WORLD,
		  &child_comm, MPI_ERRCODES_IGNORE );

  MPI_Comm_spawn( "./ERAM", argv+1,
                  num_processes_to_spawn, MPI_INFO_NULL,
                  0, MPI_COMM_WORLD,
                  &child_comm2, MPI_ERRCODES_IGNORE );
  MPI_Comm intra_comm, intra_comm2;
  MPI_Intercomm_merge( child_comm, 0, &intra_comm );
  MPI_Intercomm_merge( child_comm2, 0, &intra_comm2 );

  int intra_rank, intra_rank2;
  MPI_Comm_rank(intra_comm, &intra_rank);
  MPI_Comm_rank(intra_comm2, &intra_rank2);

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

  MPI_Comm_free(&child_comm);
  MPI_Comm_free(&intra_comm);
  MPI_Comm_free(&intra_comm2);
  MPI_Finalize();

  return 0;
}
