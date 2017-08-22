#include "petscvec.h"
#include "petsc.h"
#include "mpi.h"

/* receive information on what is the type of data that'll be received further */
int mpi_com_type_recv(com_lsa * com, int * type);
/* send to receivers what will be the type of data that'll be send */
int mpi_com_type_send(com_lsa * com, int * type);
