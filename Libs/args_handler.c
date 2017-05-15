/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "args_handler.h"

/* function used to handle component arguments this is where we choose to precond with LS-Arnoldi */
int argsHandleComponents(int argc, char ** argv, int * lsa_gmres, int * lsa_arnoldi){
	int i,rank;
	int cmpt1=0,cmpt2=0;
	
	*lsa_gmres=0;
	*lsa_arnoldi=0;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	for(i=1;i<argc;i++){
		/* do we have many components */
		if(strcmp(argv[i],"-lsa_gmres")==0){
			cmpt1=1;
			if(i+1<=argc-1){
				i++;
				*lsa_gmres=atoi(argv[i]);
			}else{
				if(rank==0)
					printf("Warning: Number of nodes for GMRES isn't set, using 1 node\n");
				*lsa_gmres=1;
			}
		}	else if(strcmp(argv[i],"-lsa_arnoldi")==0){
			cmpt2=1;
			if(i+1<=argc-1){
				i++;
				*lsa_arnoldi=atoi(argv[i]);
			}else{
				if(rank==0)
					printf("Warning: Number of nodes for ARNOLDI isn't set, using 1 node\n");
				*lsa_arnoldi=1;
			}
		}
			
	}

	if(cmpt1 && !cmpt2){
		*lsa_arnoldi=1;
		cmpt2=1;
		if(rank==0)
			printf("Warning: Number of nodes for ARNOLDI isn't set, using 1 node\n");
	}else if (!cmpt1 && cmpt2){
		*lsa_gmres=1;
		cmpt1=1;
		if(rank==0)
			printf("Warning: Number of nodes for GMRES isn't set, using 1 node\n");
	}

	if(cmpt1 && cmpt2)
		return 0;
	else 
		return 1;
}
