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

#include "triep.h"

/* 
! -----------------------------------------------------------------------
! Classes the eigenvalues in ascending order of the real parts
! Finds the rank of the last one having a real part negative
! -----------------------------------------------------------------------
 */
PetscErrorCode tri(PetscScalar * vp, PetscInt nValues, PetscInt * ch_signe){
	PetscInt i;
	/* 
	! tri par rapport aux parties reelles
	 */
	 
	for( i = 0; i < nValues-1 ; ++i)
	{
		sort(vp,0,nValues);
		
 
 		*ch_signe=0;
 	
		while(PetscRealPart(PetscRealPart(vp[*ch_signe]))<0.0){
			(*ch_signe)++;
			if(*ch_signe>nValues){
				break;
			}
		}
	}	
	
	
	return 0;
}

/* 
! Eliminates the eigenvalues whose imaginary part is negative
 */
PetscErrorCode epurer(PetscScalar * vp, PetscInt * nValues){
	PetscInt i,nKeep;
	PetscReal tmp;
	PetscScalar * eigen_keep;
	
	PetscMalloc((*nValues)*sizeof(PetscScalar),&eigen_keep);
	
	/* keep values with non negative imaginary part and non nul real part*/
	nKeep=0;
		
	for(i=0;i<*nValues;i++){
		/* abs() function does not appear as reliable, just a tweak in case of*/

///////////////

		if(PetscRealPart(vp[i])<0.0) tmp=-PetscRealPart(vp[i]);
		else tmp=PetscRealPart(vp[i]);

//////////////
		
/////		tmp=PetscRealPart(vp[i]);

		#ifdef DEBUG
		printf("$} Triep epure %e %e (%e)\n",PetscImaginaryPart(vp[i]),tmp,PetscRealPart(vp[i]));
		#endif
		if(PetscImaginaryPart(vp[i])>=0.0 && tmp!=0.0 && tmp>1.e-06){
			eigen_keep[nKeep]=vp[i];
			nKeep++;
		}		
	}
	
	for(i=0;i<nKeep;i++)
		vp[i]=eigen_keep[i];
	*nValues=nKeep;
	
	
	
	PetscFree(eigen_keep);
	
	return 0;
}

PetscErrorCode keepPositif(PetscScalar * vp, PetscInt * nValues){
  PetscInt i,nKeep;
  PetscReal tmp;
  PetscScalar * eigen_keep;

  PetscMalloc((*nValues)*sizeof(PetscScalar),&eigen_keep);

  /* keep values with non negative imaginary part and non nul real part*/
  nKeep=0;

  for(i=0;i<*nValues;i++){
    /* abs() function does not appear as reliable, just a tweak in case of*/

    ///////////////
    /*
    if(PetscRealPart(vp[i])<0.0) tmp=-PetscRealPart(vp[i]);
    else tmp=PetscRealPart(vp[i]);

    //////////////

    /////           tmp=PetscRealPart(vp[i]);

                #ifdef DEBUG
    printf("$} Triep epure %e %e (%e)\n",PetscImaginaryPart(vp[i]),tmp,PetscRealPart(vp[i]));
                #endif
    */

    if(PetscRealPart(vp[i])>0.0){
      eigen_keep[nKeep]=vp[i];
      nKeep++;
    }

  }

  for(i=0;i<nKeep;i++)
    vp[i]=eigen_keep[i];

  *nValues=nKeep;



  PetscFree(eigen_keep);

  return 0;
}

PetscErrorCode keepNegatif(PetscScalar * vp, PetscInt * nValues){
  PetscInt i,nKeep;
  PetscReal tmp;
  PetscScalar * eigen_keep;

  PetscMalloc((*nValues)*sizeof(PetscScalar),&eigen_keep);

  /* keep values with non negative imaginary part and non nul real part*/
  nKeep=0;

  for(i=0;i<*nValues;i++){
    /* abs() function does not appear as reliable, just a tweak in case of*/

    ///////////////
    /*
    if(PetscRealPart(vp[i])<0.0) tmp=-PetscRealPart(vp[i]);
    else tmp=PetscRealPart(vp[i]);

    //////////////

    /////           tmp=PetscRealPart(vp[i]);

                #ifdef DEBUG
    printf("$} Triep epure %e %e (%e)\n",PetscImaginaryPart(vp[i]),tmp,PetscRealPart(vp[i]));
                #endif
    */

    if(PetscRealPart(vp[i])<0.0){
      eigen_keep[nKeep]=vp[i];
      nKeep++;
    }

  }

  for(i=0;i<nKeep;i++)
    vp[i]=eigen_keep[i];

  *nValues=nKeep;



  PetscFree(eigen_keep);

  return 0;
}




void swap(PetscScalar *a, PetscScalar *b)
{
  PetscScalar t=*a; *a=*b; *b=t;
}


void sort(PetscScalar *arr, int beg, int end)
{
  if (end > beg + 1)
  {
    PetscScalar piv = arr[beg];
    int l ;
		int r;
    l=beg + 1;
     r = end;
    while (l < r)
    {
      if (PetscRealPart(arr[l]) <= PetscRealPart(piv))
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    sort(arr, beg, l);
    sort(arr, r, end);
  }
}










