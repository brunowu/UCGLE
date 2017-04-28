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

#include "gmres_precond.h"

static int latency_count=0;

PetscErrorCode GmresLSAPrecond(com_lsa * com, KSP ksp)
{
  KSP_FGMRES     *fgmres = (KSP_FGMRES *)(ksp->data);
  Mat            Amat,Pmat;
  PetscErrorCode ierr;
  Vec r0_tmp,w0_tmp,w1_tmp,w_1_tmp,r1_tmp,sol_tmp,vec_tmp, x_tmp;
  PetscScalar data_tmp[EIGEN_ALL*2*3];
  PetscInt size_data,ls_power,latency,hang,timing;
  PetscInt nols;
  PetscBool flag;
  int i,j, size;

  PetscScalar alpha;
  PetscScalar * eta,*delta,*beta;
  PetscReal norm;

  char  ls_load_path[PETSC_MAX_PATH_LEN];
  PetscBool ls_load, ls_load_any;

  ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_load",ls_load_path,PETSC_MAX_PATH_LEN,&ls_load);CHKERRQ(ierr);
  ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_load_any",&ls_load_any);CHKERRQ(ierr);

  if(ls_load&&ls_load_any){
	  ls_load_any=PETSC_FALSE;
	  ls_load=PETSC_TRUE;
  }

  PetscInt recv_count=0, k, tmp_size;
  PetscScalar tmp[EIGEN_ALL*2*3];

  PetscFunctionBegin;

  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_power",&ls_power,&flag);CHKERRQ(ierr);
  if(!flag) ls_power=LS_POWER;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_latency",&latency,&flag);CHKERRQ(ierr);
  if(!flag) latency=LS_LATENCY;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_m_hang",&hang,&flag);CHKERRQ(ierr);
  if(!flag) hang=ksp->max_it;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_nopc",&nols,&flag);CHKERRQ(ierr);
  if(!flag) nols=1;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_timing",&timing,&flag);CHKERRQ(ierr);
  if(!flag) timing=1000000;

  latency_count++;

  if((latency_count%latency==0  && ksp->its>0) || ksp->its==2){
    usleep(timing);
  }else return 1;

  if(nols == 0 ){
    ls_load = PETSC_FALSE;
    ls_load_any = PETSC_FALSE;
  }

  if(ls_load^=ls_load_any){
    if(!mpi_lsa_com_array_recv(com,&size_data,data_tmp)){
        PetscPrintf(PETSC_COMM_WORLD, "GMRES has recived data from LS\n");
	for(k=0;k<EIGEN_ALL*2*3;k++){
	      tmp[k] = data_tmp[k];
	    }
	    
//	PetscPrintf(PETSC_COMM_WORLD,"\n\nread\n\n");
  	PetscPrintf(PETSC_COMM_WORLD,"\n\n#}GMRES Component recived data of size %d\n\n", (PetscInt)tmp[0]);
	    tmp_size = ((PetscInt)data_tmp[0])*3+2;
	    if((((PetscInt)data_tmp[0])*3+2)!=size_data){
	      PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond data is not consistent : size mpi %d size data*3+2 %d %e data\n",
		    size_data,(((PetscInt)data_tmp[0])*3+2),((PetscReal)data_tmp[0]));
      }
    }
    else {
      for(k=0;k<EIGEN_ALL*3*2;k++){
	       data_tmp[k] = tmp[k];
	    }
	    size_data = tmp_size;
    //    PetscPrintf(PETSC_COMM_WORLD,"\n\nnot read\n\n");
    }
     //PetscPrintf(PETSC_COMM_WORLD,"\n\nnot read\n\n");
  }
  else if(nols==0||mpi_lsa_com_array_recv(com,&size_data,data_tmp)){
	   return 1;
  }

 
  //PetscPrintf(PETSC_COMM_WORLD,"\n\ntest_data=%d\n\n", size_data);
  
  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&eta);
  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&beta);
  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&delta);

  ierr=PetscMemcpy(&alpha,&data_tmp[1],1*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(eta,&data_tmp[2],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(beta,&data_tmp[2+((PetscInt)data_tmp[0])],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(delta,&data_tmp[2+2*((PetscInt)data_tmp[0])],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);

  /*get operators and vector data*/
  ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&x_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&r0_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w_1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w0_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&r1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&sol_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&vec_tmp);CHKERRQ(ierr);

  ierr=VecSet(x_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(r0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(r1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(vec_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecCopy(ksp->vec_sol,sol_tmp);CHKERRQ(ierr);
  /* proceed to computation*/

  for(j=0;j<ls_power;j++){
    /* r0 = b-Ax*/
    /* put A*x into VEC_TEMP */
    ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
    /* now put residual (-A*x + f) into vec_vv(0) */
    ierr = VecWAXPY(r0_tmp,-1.0,vec_tmp,ksp->vec_rhs);CHKERRQ(ierr);
    /* r0 = w0*/
    ierr=VecCopy(r0_tmp,w0_tmp);CHKERRQ(ierr);
    ierr=VecCopy(w0_tmp,x_tmp);CHKERRQ(ierr);
    ierr=VecScale(x_tmp,eta[0]);CHKERRQ(ierr);
    ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
    /* depending ton the ls polynom size (PetscInt)data_tmp[0] */
    for(i=0;i<(PetscInt)data_tmp[0]-1;i++){
      /* w1=-alpha*w0 - delta[i]*w_1 ((  y = alpha x + delta y. )) (Vec y,PetscScalar alpha,PetscScalar beta,Vec x)*/
      ierr=VecCopy(w_1_tmp,w1_tmp);CHKERRQ(ierr);
      ierr=VecAXPBY(w1_tmp,-alpha,-(PetscScalar)delta[i],w0_tmp);CHKERRQ(ierr);
      /* w1 = w1 - A*w0 */
      ierr = MatMult(Amat,w0_tmp,vec_tmp);CHKERRQ(ierr);
      /* y = alpha x + y.  VecAXPY(Vec y,PetscScalar alpha,Vec x)*/
      ierr = VecAXPY(w1_tmp,1.0,vec_tmp);CHKERRQ(ierr);
      /* w1 = w1/beta[i] */
      ierr=VecScale(w1_tmp,1/beta[i]);CHKERRQ(ierr);
      /* w_1 = w0 */
      ierr=VecCopy(w0_tmp,w_1_tmp);CHKERRQ(ierr);
      /* w0 = w1 */
      ierr=VecCopy(w1_tmp,w0_tmp);CHKERRQ(ierr);
      /* x = x + (w0 * eta[i] ) */
      ierr=VecAXPY(x_tmp,eta[i+1],w0_tmp);CHKERRQ(ierr);
    }
    /* x1= x1+x*/
    ierr=VecAXPY(sol_tmp,1.0,x_tmp);CHKERRQ(ierr);
    /* put A*x into VEC_TEMP */
    ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
    /* now put residual (-A*x + f) into vec_vv(0) */
    ierr = VecWAXPY(r1_tmp,-1.0,vec_tmp,ksp->vec_rhs);CHKERRQ(ierr);
    /* compute norm and see if it's below epsilon */
    VecNorm(r1_tmp,NORM_2,&norm);
  }

  if(nols!=0) ierr=VecCopy(sol_tmp,ksp->vec_sol);CHKERRQ(ierr);
  ierr=PetscFree(eta);CHKERRQ(ierr);
  ierr=PetscFree(beta);CHKERRQ(ierr);
  ierr=PetscFree(delta);CHKERRQ(ierr);
  return 0;
}
