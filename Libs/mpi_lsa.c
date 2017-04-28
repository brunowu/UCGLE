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

#include "mpi_lsa.h"

/* init lsa mpi layer if needed */
int mpi_lsa_init(int argc, char ** argv, com_lsa * com){
	/* number of nodes for each lsa component */
	int lsa_gmres;
	int lsa_arnoldi;
	int lsa_ls=1;
	int lsa_father=1;
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	/* before initializing anything we gonna check if we need to use lsa */
	if(argsHandleComponents(argc, argv, &lsa_gmres, &lsa_arnoldi)){
		if(rank==0)
			printf("NOTE : Using classical GMRES, no LSA mode.\n");
		return 1;
	}
	if(rank==0)
		printf("NOTE : Using LSA mode.\n");

	/* get some general information */
	com->com_world=MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD,&(com->rank_world));
	MPI_Comm_size(MPI_COMM_WORLD,&(com->size_world));

	/* put the groups sizes into lsa struct
	choice of ordering is arbitrarly made */
	com->size.com[0]=lsa_gmres;	//ici lsa_gmres est non-initialisé
	com->size.com[1]=lsa_father;
	com->size.com[2]=lsa_arnoldi;	//ici lsa_arnoldi est non-initialisé
	com->size.com[3]=lsa_ls;
	com->nbr_array_sended = 0;
	if(com->rank_world==0)printf("]> Creating Groups\n");
	/* now we need to generate the mpi communicators for each group */
	mpi_lsa_create_groups(com);
	/* and create communicators between the groups */
	mpi_lsa_create_intercoms(com);

	if(com->rank_world==0)printf("]> MPI communications ready to be used\n");
	return 0;
}
/* create lsa groups and communicators */
int mpi_lsa_create_groups(com_lsa * com){
	int i,j,tmp_int;

	/* get the color */
	/* get the id of the first process for a group
	 Here we just retrieve the rank of the master (first of the group) for each group (rank in MPI_COMM_WORLD) */
	for(i=0;i<4;i++){
		tmp_int=0;
		for(j=0;j<i;j++){
			tmp_int+=com->size.com[j];
		}
		com->master.com[i]=tmp_int;
	}
	tmp_int=0;
	for(i=0;i<4;i++){
		if(com->rank_world>=com->master.com[i]){
			com->color_group=i;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/* create a communicator for each group of subprocess */
	/* we create as many groups as there are distinct values for color
	 	 the rank of the process in the group depends on its rank in comm_wold "com->rank_world"
		 the group to wich belongs the process is saved in com->com_group */
	MPI_Comm_split(MPI_COMM_WORLD,com->color_group,com->rank_world,&com->com_group);
	MPI_Comm_rank(com->com_group,&com->rank_group);
	MPI_Comm_size(com->com_group,&tmp_int);	// possibible to collect the size of com_group

	return 0;
}

/* create intercommunicators between intracoms */
/* create a loop between all the elements types */
int mpi_lsa_create_intercoms(com_lsa * com){
	int prev, next,flag;
	int prev_size,next_size,size;
	/* create first connection between intracommunicators thanks to an intercommunicator */
	/* one way */
	MPI_Barrier(MPI_COMM_WORLD);
	if(com->rank_world==0)printf("]> Creating intercommunicators\n-One Way :\n");
	MPI_Barrier(MPI_COMM_WORLD);
	printf("\t *> %d -> %d ",com->rank_world,com->master.com[4-((com->color_group)+1)]);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Intercomm_create(com->com_group,0,
					MPI_COMM_WORLD,com->master.com[4-((com->color_group)+1)],
					com->rank_group,
					&(com->inter.com[4-((com->color_group)+1)]));

	MPI_Barrier(MPI_COMM_WORLD);
	if(com->rank_world==0)printf("\n]> The Other : \n");
	MPI_Barrier(MPI_COMM_WORLD);
	printf("\t *> %d -> %d ",(com->color_group),(4-((com->color_group)-1)%4)%4);
	MPI_Barrier(MPI_COMM_WORLD);
	if(com->rank_world==0)printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);

 	/* the other */
	MPI_Intercomm_create(com->com_group,0,
					com->com_world,com->master.com[(4-((com->color_group)-1)%4)%4],
	 				com->rank_group,
	 				&(com->inter.com[(4-((com->color_group)-1)%4)%4]));

	if((4-((com->color_group)-1)%4)%4>com->color_group){
		next=(4-(com->color_group-1)%4)%4;
		prev=4-((com->color_group)+1);
	} else {
		prev=(4-(com->color_group-1)%4)%4;
		next=4-((com->color_group)+1);
	}

	/* set the in and out communicators */
	com->out_com=com->inter.com[next];
	MPI_Comm_test_inter(com->inter.com[next],&flag);
		if(!flag){
			mpi_lsa_print("\n\n\n\nproblem with inter.[next]\n\n\n\n\n", com);
		}
	com->in_com=com->inter.com[prev];
	MPI_Comm_test_inter(com->inter.com[prev],&flag);
	if(!flag){
		mpi_lsa_print("\n\n\n\n\nproblem with inter.[prev]\n\n\n\n", com);
	}

	MPI_Comm_remote_size(com->out_com,&next_size);
	MPI_Comm_remote_size(com->in_com,&prev_size);
	MPI_Comm_size(com->com_group,&size);
	if(com->rank_world==0) printf("]> In and Out communicators : \n");
	MPI_Barrier(MPI_COMM_WORLD);

	if(com->color_group==0) 		 printf("GMRES :   ");
	else if(com->color_group==1) printf("MAIN :    ");
	else if(com->color_group==2) printf("ARNOLDI : ");
	else if(com->color_group==3) printf("LS :      ");

	printf("%d: %d (%d) -> %d (%d) -> %d (%d)\n",com->rank_world,com->master.com[prev],prev_size,com->color_group,size,com->master.com[next],next_size);


	return 0;
}



int mpi_lsa_print(char * s,com_lsa * com){

	if(com->rank_world==0) printf("\n");
	return 0;
}
