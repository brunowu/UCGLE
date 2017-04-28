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

 #include "lsqr.h"

PetscErrorCode LSQR(com_lsa * com, int * vector_size){
	/* variables */
	PetscInt end,cumul,eigen_received,eigen_total,eigen_max;
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	int i,info,type=0;
	PetscBool flag,data_load,data_export,continuous_export,data_load_any;
	PetscScalar * data,*eigen_cumul,*eigen_tri,*d,*c;
	PetscReal a_ell,c_ell,d_ell,d_reel;
	int data_size;
	PetscInt chsign;
	PetscInt mu1,mu2,mu,result_array_size;
	PetscScalar alpha,*eta,*beta,*delta,scalar_tmp;
	PetscInt ls_eigen_min, ls_eigen; // use default values
	PetscErrorCode ierr;
 	PetscScalar * result_array,*data_buffer;/*[EIGEN_ALL*3+2]*/
	Vec * v;
	sprintf(load_path,"./lsqr.bin");
	sprintf(export_path,"./lsqr.bin");
	/* check if there is arguments for ls */
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen_min",&ls_eigen_min,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen_min=EIGEN_MIN;
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen",&ls_eigen,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen=EIGEN_ALL;
	/* check the number of eigenvalues that one will receive from arnoldi */
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_k_param",&eigen_max,&flag);CHKERRQ(ierr);
	if(!flag)eigen_max=ls_eigen;

	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_load_any",&data_load_any);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_cexport",&continuous_export);CHKERRQ(ierr);

	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&eigen_tri);
	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&eigen_cumul);
	ierr=PetscMalloc(((*vector_size)+1)*sizeof(PetscScalar),&d);
	ierr=PetscMalloc(((*vector_size)+1)*sizeof(PetscScalar),&c);
	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&data);
	ierr=PetscMalloc(eigen_max*sizeof(PetscScalar),&data_buffer);

	/* data that will be sended to GMRES for it's preconditionning step */
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&eta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&beta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&delta);
	ierr=PetscMalloc((eigen_max*3+2)*sizeof(PetscScalar),&result_array);
	result_array_size=2+3*eigen_max;

	for(i=0;i<(*vector_size);i++){
	  eigen_tri[i]=(PetscScalar)0.0;
	  eigen_cumul[i]=(PetscScalar)0.0;
	}
	for(i=0;i<(*vector_size)+1;i++){
	  d[i]=(PetscScalar)0.0;
	  c[i]=(PetscScalar)0.0;
	}
	cumul=0;
	eigen_received=0;
	eigen_total=0;
	end=0;
	eigen_total=0;
	ls_eigen=0;
	while(!end){
		if(!mpi_lsa_com_type_recv(com,&type)){
		  if(type==666){
		    end=1;
		    break;
		  }else if (type == 911){
		  	mpi_lsa_com_type_send(com,&type);
		  }
		}
		/*in any case clear data array*/
		for(i=0;i<eigen_max;i++)
		  data[i]=(PetscScalar)0.0+PETSC_i*(PetscScalar)0.0;
		/* if received something */
		if(!mpi_lsa_com_array_recv(com, &data_size,data) || data_load || data_load_any){
			/* we received data or load it depending on the flags (for first step only*/
			if(data_load&&data_load_any){
			  data_load_any=PETSC_FALSE;
			  data_load=PETSC_TRUE;
			}
			if(!(data_load^=data_load_any)){
				/* first we gonna remove some non-needed values */
				epurer(data,&data_size);
				/* add them to the accumulated eigenvalues */
				/* if full renew full eigenvalues */
				if(eigen_total+data_size>*vector_size) eigen_total=0;
				/* select eigenvalues */
				for(i=0;i<data_size;i++){
					eigen_cumul[eigen_total+i]=data[i];
				}
				eigen_total+=data_size;
				if(cumul<eigen_total) cumul=eigen_total;
				for(i=0;i<cumul;i++){
					eigen_tri[i]=eigen_cumul[i];
				}
			} else {
				ierr=readBinaryScalarArray(load_path,&cumul, eigen_tri);CHKERRQ(ierr);
				PetscPrintf(PETSC_COMM_WORLD, "LS has read the local file\n");
				data_load=PETSC_FALSE;
				data_load_any=PETSC_FALSE;
				data_size=cumul;
			}
			eigen_received+=data_size;
			/* if we didn't received enough eigenvalues */
			if(eigen_received<ls_eigen_min && data) continue;
			else {
				eigen_received=0;
				tri(eigen_tri,cumul,&chsign);
				mu1=0;
				mu2=0;
				/* convex hull computation */
				if(chsign>0){
					//keepPositif(eigen_tri,&cumul);
					convhull(eigen_tri, c, d, chsign, &mu1, 0, 0);
					printf("@} LSQR convhul negatif chsigne %d cumul %d mu1 %d\n",chsign,cumul,mu1);

				}
				if(chsign<cumul){
					convhull(eigen_tri, c, d, cumul-chsign, &mu2, chsign, mu1);
					printf("@} LSQR convhul positif chsigne %d cumul %d mu1 %d mu2 %d\n",chsign,cumul,mu1,mu2);

				}
				mu=mu1+mu2;
				/* Ellipse computation */
				ierr=ellipse(c,  d, mu+1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info);CHKERRQ(ierr);
				if(fabs(d_ell)<epsilon()) d_ell = 1.;
				if(fabs(a_ell)<epsilon())
					ls_eigen=0;
				else {
					LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta,
					  delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max);
				}
			}
		}

		if(ls_eigen>1){
			/* place the computed results inside the array */
			scalar_tmp=(PetscScalar)ls_eigen;
			ierr=PetscMemcpy(&result_array[0],&scalar_tmp,1*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[1],&alpha,1*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2],eta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+ls_eigen],beta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+(2*ls_eigen)],delta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);
			if(continuous_export){
			  ierr=writeBinaryScalarArray(data_export?export_path:"./lsqr.bin", cumul, eigen_tri);
			}
			/* and send it */
			result_array_size=2+3*ls_eigen;
			mpi_lsa_com_array_send(com, &result_array_size,result_array);
                        PetscPrintf(PETSC_COMM_WORLD, "LS has sent the parameters\n");
		}
		if(ls_eigen>1){
		  ls_eigen=0;
		}
	}
	if(data_export){
		ierr=writeBinaryScalarArray(export_path, cumul, eigen_tri);
	}
	/* Free the arrays */
	PetscFree(eigen_tri);
	PetscFree(eigen_cumul);
	PetscFree(d);
	PetscFree(c);
	PetscFree(data);
	PetscFree(eta);
	PetscFree(beta);
	PetscFree(delta);
 	ierr=PetscFree(result_array);CHKERRQ(ierr);
	return 0;
}
