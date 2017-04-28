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

#ifndef MPI_LSA_H_
#define MPI_LSA_H_

#include "unistd.h"
#include "stdio.h"
#include "mpi.h"
#include "args_handler.h"
#include "petsc.h"

/* some constants */
#define GMRES_GR 0
#define FATHER_GR 1
#define ARNOLDI_GR 2
#define LS_GR 3


typedef struct _com_lsa{
	int com_world;							// MPI_COMM_WORLD

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	}group;

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} inter; 								// intercommucators: it contains the intercomm with the leader of the peer group

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} size;									// groups size: it indicates the size of each group

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} master;								// group-masters: it indicates the rank_world of the master of each group

	int size_world;						// size of MPI_COMM_WORLD
	int com_group;							// group's communicator
	int rank_group;						//rank in the group
	int color_group;						// the color of the goup to which belongs the process

	int rank_world;						// rank in world
	int recv_flag;							// indicates if there a msg to be received by the group( the master get the information and distribut it)
												// if true then the master will receive it and the others will wait for a scatter.

	int * in_size;							//
	int * out_size;						//

//	requests array for outcomming messaging : their size depends on out_number (see bellow)
	int * vec_requests;					// this for sendig Petsc Vec
	int * type_requests;					// this for sending special messages like exit message
	int * array_requests;				// this for sending arrays

	int vec_requests_nb;
	int type_requests_nb;
	int array_requests_nb;

	int * vec_in_disp;					// table of displacement inside the incomming vec
	int * vec_out_disp;					// array of displacement inside the outcomming vec
	int * vec_color_disp;

	int in_number;							// size of in_comm: means the number of process in the incomming comm group
	int out_number;						// size of out_com: means the number of process in the outcomming comm group

	int out_sended;						// counter of sended msg
	int in_received;						// counter of received msg

	PetscScalar * out_sended_buffer;
	PetscScalar * in_received_buffer;

	int nbr_array_sended;
	PetscScalar * array_out_sended_buffer;
	PetscScalar * array_in_received_buffer;

	int out_vec_sended;					// counter of sended vecs ( including those sent and not received by the remot process)
	int array_out_sended;				// counter of sended arrays ( including those sent and not received by the remot process)

	int in_com;								// communicator for incomming messaging
	int out_com;  							// communicator for outputing messaging

	int * vec_color_sizes;				//Vector sizes within the group: means the size of the vector for each process in the color_group
	int * vec_out_sizes;					//Vector sizes within the out group: means the size of the vector for each process in the out_group
	int * vec_in_sizes;					//Vector sizes within the in group: means the size of the vector for each process in the in_group

} com_lsa;

int mpi_lsa_init(int argc, char ** argv, com_lsa * com);

int mpi_lsa_create_groups(com_lsa * com);

int mpi_lsa_create_intercoms(com_lsa * com);

int mpi_lsa_print(char * s,com_lsa * com);

#endif
