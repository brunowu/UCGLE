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

#include "ellipse.h"

int ellipse3(PetscScalar xy1, PetscScalar xy2, PetscScalar xy3, PetscReal * a2, PetscReal * b2, PetscReal * c, PetscInt * info){
	
	PetscReal d2,z,r;
			
	r = (PetscRealPart(xy2)-PetscRealPart(xy1))*(PetscRealPart(xy3)-
			 PetscRealPart(xy2))*(PetscRealPart(xy3)-PetscRealPart(xy1));
	
	if(r>1.e-14){
		z=pow(PetscImaginaryPart(xy1),2)*(PetscRealPart(xy2)-PetscRealPart(xy3))+
			pow(PetscImaginaryPart(xy2),2)*(PetscRealPart(xy3)-PetscRealPart(xy1))+
			pow(PetscImaginaryPart(xy3),2)*(PetscRealPart(xy1)-PetscRealPart(xy2));

		if(fabs(z)<1.e-14) {
		  return 1;
		}
		  
		*c=(pow(PetscImaginaryPart(xy1),2)*(pow(PetscRealPart(xy2),2)-pow(PetscRealPart(xy3),2))+
			  pow(PetscImaginaryPart(xy2),2)*(pow(PetscRealPart(xy3),2)-pow(PetscRealPart(xy1),2))+
				pow(PetscImaginaryPart(xy3),2)*(pow(PetscRealPart(xy1),2)-pow(PetscRealPart(xy2),2)))/(2*z);
			
					
					
		*a2=pow(*c,2)-(pow(PetscImaginaryPart(xy1),2)*PetscRealPart(xy2)*PetscRealPart(xy3)*(PetscRealPart(xy2)-PetscRealPart(xy3))
			     +pow(PetscImaginaryPart(xy2),2)*PetscRealPart(xy1)*PetscRealPart(xy3)*(PetscRealPart(xy3)-PetscRealPart(xy1))
			     +pow(PetscImaginaryPart(xy3),2)*PetscRealPart(xy1)*PetscRealPart(xy2)*(PetscRealPart(xy1)-PetscRealPart(xy2)))/z;
		
		
		d2=(*a2)*(1.0 - z/r);

		(*b2)=(*a2)-d2;

	} else {
		return 1;
	}
	
	*info=1;
	if(*a2<0.0 || *b2<0.0) *info=0;																								
	return 0;												
}


int ellipse2(PetscScalar xy1, PetscScalar xy2, PetscReal * a2, PetscReal * b2, PetscReal * c, PetscInt * info){
	PetscReal a,b,s,t,z,q,d2,sign;												
	
	a=(PetscRealPart(xy2)-PetscRealPart(xy1))/2;
	b=(PetscRealPart(xy2)+PetscRealPart(xy1))/2;
	s=(PetscImaginaryPart(xy2)-PetscImaginaryPart(xy1))/2;
	t=(PetscImaginaryPart(xy2)+PetscImaginaryPart(xy1))/2;
	
	if(t<=1.e-15 || fabs(a)<=1.e-15)	{
	  return 1;
	}
	
	if(fabs(s)>=1.e-15){
		q=(s/t+t/s)/2;
		sign=1.0;

		if(a*s<=0) sign=-1.0;
		
		z=((-q)+sign*sqrt(pow(q,2)+3))*(a/3);
		if(fabs(z)<1.e-15){ 
		  return 1;
		}
		*a2=(z + a*t/s) * (z + a*s/t);
		d2=(*a2)*(z - s*t/a)/z;
		*b2=(*a2)-d2;
		*c=b+z;
	} else {
		*a2=2*pow(a,2);
		*b2=2*pow(*c,2);
		d2=(*a2)-(*b2);
		z=0.0;
		*c=b+z;		
	}			
		
	*info=1;
	if((*a2)<0.0 || (*b2)<0.0) *info=0;
	

	return 0;
}


int test(PetscScalar * hk, PetscInt n, PetscReal a2, PetscReal b2, PetscReal c, PetscInt * info){
	PetscInt i;
	PetscReal d;
	
	*info=1;
	for(i=0;i<n;i++){
		d=b2*pow(PetscRealPart(hk[i])-c,2)+a2*pow(PetscImaginaryPart(hk[i]),2)-a2*b2;
		if(d>=1.e-6){
			*info=0;
			return 1;
		}
	}
	
	return 0;
}

PetscErrorCode ellipse(PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt mu, PetscReal * co, PetscReal * ao2, PetscReal * do2, PetscReal * dr, PetscInt * info){
	PetscInt i,j,k,info1,info2;
	PetscScalar * hk;
	PetscErrorCode ierr;
	PetscReal bo2,a2,b2,cc,aire;
	
	#ifdef DEBUG
		PetscPrintf(PETSC_COMM_WORLD,"$} Ellipse Allocating work memory %d\n",n);
	#endif
	ierr=PetscMalloc(n*sizeof(PetscScalar),&hk);CHKERRQ(ierr);
	
	#ifdef DEBUG
		PetscPrintf(PETSC_COMM_WORLD,"$} Ellipse Computing Edges\n");
	#endif

	for(i=0;i<n-1;i++){
		hk[i]=c[i]-d[i];
	}
	
	if(mu!=0 && mu<n){
		hk[n-1]=c[mu]+d[mu];
	}
	
	hk[n-1]=c[n-2]+d[n-2];
	
	
	#ifdef DEBUG
		PetscPrintf(PETSC_COMM_WORLD,"$} Ellipse Research two points optimal ellipse\n");
	#endif
	
	/* find two points optimal ellipse */
	info1=0;
	info2=0;
	
	*ao2=0.0;
	bo2=0.0;
	*co=0.0;
	
	aire=1.e16;
	
	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			info1=0;
			ellipse2(hk[i],hk[j],&a2,&b2,&cc,&info1);
			if(info1==0) continue;
			test(hk,n,a2,b2,cc,&info1);
			
			if(info1==1){
				if(a2*b2<aire){
					info2=1;
					(*ao2)=a2;
					bo2=b2;
					(*co)=cc;
					aire=(*ao2)*bo2;		
				}
			}
		}
	}
	
	#ifdef DEBUG
		PetscPrintf(PETSC_COMM_WORLD,"$} Ellipse Research three points optimal ellipse\n");
	#endif
	
	/* find three points optimal ellipse */
	
	if(info2!=1){
		aire=1.e16;
		
		for(i=0;i<n-2;i++){
			for(j=i+1;j<n-1;j++){
				for(k=j+1;k<n;k++){
					*info=0;					
					ellipse3(hk[i],hk[j],hk[k],&a2,&b2,&cc,info);


					if(*info==0) continue;
					test(hk,n,a2,b2,cc,info);
					if(*info==1){
						if(a2*b2<aire){
							*ao2=a2;
							bo2=b2;
							*co=cc;
							aire=(*ao2)*b2;							
						}
					}
				}
			}
		}
	}
	
	*do2=(*ao2)-bo2;
	*dr=(*ao2)>=bo2;
	*ao2=sqrt(*ao2);
	
	ierr=PetscFree(hk);CHKERRQ(ierr);
	
	return 0;
}









