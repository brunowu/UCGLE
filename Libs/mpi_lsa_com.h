/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#ifndef MPI_LSA_COM_H_
#define MPI_LSA_COM_H_


#include "petscvec.h"
#include "petsc.h"
#include "mpi.h"
#include "mpi_lsa.h"

#ifndef T_VEC
#define T_VEC 1
#endif

#ifndef T_TAB
#define T_TAB 2
#endif

#ifndef T_CTRL
#define T_CTRL 0
#endif


/* exchange the vector repartition for send and receive purpose */
int mpi_lsa_com_vecsize_init(com_lsa * com, Vec * v);



/* send a vector */
int mpi_lsa_com_vec_send(com_lsa * com, Vec * v);
/* receive a vector */
int mpi_lsa_com_vec_recv(com_lsa * com, Vec * v);
/* validate if a vector has entirely been received */
int mpi_lsa_com_vec_recv_validate(com_lsa * com, Vec * v, int size);



/* receive information on what is the type of data that'll be received further */
int mpi_lsa_com_type_recv(com_lsa * com, int * type);
/* send to receivers what will be the type of data that'll be send */
int mpi_lsa_com_type_send(com_lsa * com, int * type);

int mpi_lsa_send_vec(com_lsa * com, Vec * v);

int mpi_lsa_receive_vec(com_lsa * com, Vec * v);

/* send an array of data of PetscScalar type */
int mpi_lsa_com_array_send(com_lsa * com, int * size, PetscScalar * data);
/* receive an array of data of PetscScalar type, data array must be allocated */
int mpi_lsa_com_array_recv(com_lsa * com, int * size, PetscScalar * data);

int setting_out_vec_sizes( com_lsa * com, Vec * v);

int mpi_lsa_com_free(com_lsa * com);

#endif
