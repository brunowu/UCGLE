/*Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "convhul.h"

int convhull(PetscScalar * ab, PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt * ne, PetscInt offset, PetscInt mu){
	PetscScalar * hk;
	PetscReal s,cs,dzero,deps;
	PetscInt m,i,j,m1;
	PetscInt * l;
	PetscErrorCode ierr;
	
	/* Offset:  Acces aux elements [1+offset:n+offset]
                     ! des tableaux a et b */
	/* mu: 	! Acces aux elements [1+mu:n+mu]
	 ! des tableaux cr, ci, dr, di */
	
	ierr=PetscMalloc((n+2)*sizeof(PetscInt),&l);CHKERRQ(ierr);
	ierr=PetscMalloc((n+2)*sizeof(PetscScalar),&hk);CHKERRQ(ierr);
	
	dzero=0.0;
	deps=1.e-16;
	
	for(i=0;i<n+2;i++){
		l[i]=1;
	}
	
	
	/* trouve l'indice de la valeur propres à plus petite partie réelle*/
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
	
	/* met la partie imaginaire de la valeur propre de plus petite partie réelle à 0 */
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





