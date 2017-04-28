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

#include "gmres.h"

PetscErrorCode launchGMRES(com_lsa * com, Vec * b, Mat * A){
	PetscErrorCode ierr;
	KSP ksp, kspp;
	Vec x, c;
	KSPConvergedReason reason;
	PetscInt its, nols, ntimes;
	int i, size;
	PetscBool flagls, flagtimes;
	double cost_time;
	clock_t start, end;

	PetscReal norm;
	VecGetSize(*b, &size);
	VecDuplicate(*b, &c);
        PetscPrintf(com->com_group,"#} GMRES Creating and setting vector x\n");

//	generate_random_seed_vector(size, -10,10, 10,&c);
	ierr = VecDuplicate(*b, &x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSet(x, 0.0); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_nopc",&nols,&flagls);CHKERRQ(ierr);
	ierr = KSPCreate(com->com_group, &ksp);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
//	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	ierr = KSPSetOperators(ksp, *A, *A);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);


//	PetscPrintf(com->com_group,"#} GMRES Creating and setting vector x\n");

	PetscOptionsGetInt(NULL,"-ntimes",&ntimes,&flagtimes);

	if(!flagtimes){
		ntimes=1;
	}

	for (i=1; i<=ntimes; i++) {
		if(i==1){
			VecCopy(*b,c);
			VecNorm(c, NORM_2,&norm);
			VecScale(c, 0.01/norm);
		}
		else{
			generate_random_seed_vector(size, -10,10, i,&c);
			VecNorm(c, NORM_2,&norm);
			VecScale(c, 0.01/norm);
		}
		start=clock();
//		ierr = KSPSolve(ksp, *b, x); CHKERRQ(ierr);

		ierr = MyKSPSolve(ksp, c, x,com); CHKERRQ(ierr);
		end=clock();
  		cost_time = (double)(end - start)/CLOCKS_PER_SEC;
		KSPGetConvergedReason(ksp,&reason);
		if (reason<0) {
			PetscPrintf(PETSC_COMM_WORLD,"\nResolution %d: Divergence in acceptable iteration steps.\n", i);
		}
		else {
			KSPGetIterationNumber(ksp,&its);
			PetscPrintf(PETSC_COMM_WORLD,"\nResolution %d: Convergence in %f seconds / %d iterations. \n", i, cost_time, its);
		}
		}
	int exit_type=666;
	mpi_lsa_com_type_send(com,&exit_type);
	PetscPrintf(com->com_group,"\n\n#}Finish the resolution\n");
  	return ierr;
}
