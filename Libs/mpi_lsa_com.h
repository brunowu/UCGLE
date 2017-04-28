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
