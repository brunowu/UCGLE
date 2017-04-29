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

#include "convhul.h"

int convhull(PetscScalar * ab, PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt * ne, PetscInt offset, PetscInt mu){
	PetscScalar * hk;
	PetscReal s,cs,dzero,deps;
	PetscInt m,i,j,m1;
	PetscInt * l;
	PetscErrorCode ierr;
	
	ierr=PetscMalloc((n+2)*sizeof(PetscInt),&l);CHKERRQ(ierr);
	ierr=PetscMalloc((n+2)*sizeof(PetscScalar),&hk);CHKERRQ(ierr);
	
	dzero=0.0;
	deps=1.e-16;
	
	for(i=0;i<n+2;i++){
		l[i]=1;
	}
	
	/* find the index of eigenvalues with most small real part*/
	s=(PetscScalar)PetscRealPart(ab[offset])+PETSC_i*(PetscScalar)0.0;
	m1=0;
	
	for(i=1;i<n;i++){
		if(PetscRealPart(ab[i+offset])<=PetscRealPart(s)){
			s=(PetscScalar)PetscRealPart(ab[offset])+PETSC_i*(PetscScalar)0.0;
	       		
///////			s=(PetscScalar)PetscRealPart(ab[i+offset])+PETSC_i*(PetscScalar)0.0;
			m1=i;
//////			PetscPrintf(PETSC_COMM_WORLD, "convex hull test points\n");
		}
	}
	
	*ne=0;
	hk[0]=s+PETSC_i*dzero; /* init complex array by s+i*0 */
	
	if(PetscImaginaryPart(ab[m1+offset])<deps){
		l[m1]=0;
	}else{
////////		m1=0;
   		ab[m1+offset]=(PetscScalar) PetscRealPart(s) + (PetscScalar)PETSC_i*dzero;
	}
	
	
	for( i = 1; i <= n; ++i){
		s=-1.0;
		m=m1;
		for( j = 0; j < n; ++j){
			if(l[j]){
				cs=sqrt(pow(PetscRealPart(ab[j+offset])-PetscRealPart(hk[i-1]),2) +
				       (pow(PetscImaginaryPart(ab[j+offset])-PetscImaginaryPart(hk[i-1]),2)));
				if(cs>=deps){
					cs=(PetscImaginaryPart(ab[j+offset])-PetscImaginaryPart(hk[i-1]))/cs;
				}else{
					continue;
				}
				if(cs-s>-deps){
					if(cs-s<deps){
						if(PetscRealPart(ab[j+offset]) > PetscRealPart(ab[m+offset]+deps)){ 
						  m=j;
						}
						if(s>=(1-deps) && PetscImaginaryPart(ab[j+offset])>PetscImaginaryPart(ab[m+offset])){
						  m=j;
						}
					}else {
						s=cs;
						m=j;
					}
				}
			}
		}
		if(m==m1){
			break;	
		}else {
			hk[i]=ab[m+offset];
			l[m]=0;
			*ne=i;
			c[i-1+mu]=(hk[i]+hk[i-1])/2;
			d[i-1+mu]=(hk[i]-hk[i-1])/2;
			
			for(j=0;j<n;j++)
				if(l[j])
					if(PetscRealPart(ab[j+offset])<=PetscRealPart(hk[i])) l[j]=0;
		}
	}
	
	if(PetscImaginaryPart(hk[*ne])!=dzero){
		c[(*ne)+mu]=PetscRealPart(hk[*ne])+PetscImaginaryPart(hk[*ne])/2*PETSC_i;
		d[(*ne)+mu]=dzero+PetscImaginaryPart(hk[*ne])/2*PETSC_i;
		*ne=(*ne)+1;
	}

	PetscFree(l);
	PetscFree(hk);
	
	return 0;
}





