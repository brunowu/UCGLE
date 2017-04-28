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

#include "gmres_cycle.h"

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESCycle"
PetscErrorCode MyKSPFGMRESCycle(PetscInt *itcount,KSP ksp)
{

  KSP_FGMRES     *fgmres = (KSP_FGMRES*)(ksp->data);
  PetscReal      res_norm;
  PetscReal      hapbnd,tt;
  PetscBool      hapend = PETSC_FALSE;  /* indicates happy breakdown ending */
  PetscErrorCode ierr;
  PetscInt       loc_it;                /* local count of # of dir. in Krylov space */
  PetscInt       max_k = fgmres->max_k; /* max # of directions Krylov space */
  Mat            Amat,Pmat;

	PetscFunctionBegin;
  /* Number of pseudo iterations since last restart is the number
     of prestart directions */
  loc_it = 0;

  /* note: (fgmres->it) is always set one less than (loc_it) It is used in
     KSPBUILDSolution_FGMRES, where it is passed to KSPFGMRESBuildSoln.
     Note that when KSPFGMRESBuildSoln is called from this function,
     (loc_it -1) is passed, so the two are equivalent */
  fgmres->it = (loc_it - 1);

  /* initial residual is in VEC_VV(0)  - compute its norm*/
  ierr = VecNorm(VEC_VV(0),NORM_2,&res_norm);CHKERRQ(ierr);

  /* first entry in right-hand-side of hessenberg system is just
     the initial residual norm */
  *RS(0) = res_norm;

  ksp->rnorm = res_norm;
  ierr       = KSPLogResidualHistory(ksp,res_norm);CHKERRQ(ierr);
  ierr       = KSPMonitor(ksp,ksp->its,res_norm);CHKERRQ(ierr);

  /* check for the convergence - maybe the current guess is good enough */
  ierr = (*ksp->converged)(ksp,ksp->its,res_norm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
  if (ksp->reason) {
    if (itcount) *itcount = 0;
    PetscFunctionReturn(0);
  }

  /* scale VEC_VV (the initial residual) */
  ierr = VecScale(VEC_VV(0),1.0/res_norm);CHKERRQ(ierr);

  /* MAIN ITERATION LOOP BEGINNING*/
  /* keep iterating until we have converged OR generated the max number
     of directions OR reached the max number of iterations for the method */
  while (!ksp->reason && loc_it < max_k && ksp->its < ksp->max_it) {
    if (loc_it) {
      ierr = KSPLogResidualHistory(ksp,res_norm);CHKERRQ(ierr);
      ierr = KSPMonitor(ksp,ksp->its,res_norm);CHKERRQ(ierr);
    }
    fgmres->it = (loc_it - 1);

    /* see if more space is needed for work vectors */
    if (fgmres->vv_allocated <= loc_it + VEC_OFFSET + 1) {
      ierr = MyKSPFGMRESGetNewVectors(ksp,loc_it+1);CHKERRQ(ierr);
      /* (loc_it+1) is passed in as number of the first vector that should
         be allocated */
    }
    /* CHANGE THE PRECONDITIONER? */
    /* ModifyPC is the callback function that can be used to
       change the PC or its attributes before its applied */
    (*fgmres->modifypc)(ksp,ksp->its,loc_it,res_norm,fgmres->modifyctx);

    /* apply PRECONDITIONER to direction vector and store with
       preconditioned vectors in prevec */
    ierr = KSP_PCApply(ksp,VEC_VV(loc_it),PREVEC(loc_it));CHKERRQ(ierr);
    ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);
    /* Multiply preconditioned vector by operator - put in VEC_VV(loc_it+1) */
    ierr = KSP_MatMult(ksp,Amat,PREVEC(loc_it),VEC_VV(1+loc_it));CHKERRQ(ierr);
    /* update hessenberg matrix and do Gram-Schmidt - new direction is in
       VEC_VV(1+loc_it)*/
    ierr = (*fgmres->orthog)(ksp,loc_it);CHKERRQ(ierr);
    /* new entry in hessenburg is the 2-norm of our new direction */
    ierr = VecNorm(VEC_VV(loc_it+1),NORM_2,&tt);CHKERRQ(ierr);
    *HH(loc_it+1,loc_it)  = tt;
    *HES(loc_it+1,loc_it) = tt;
    /* Happy Breakdown Check */
    hapbnd = PetscAbsScalar((tt) / *RS(loc_it));
    /* RS(loc_it) contains the res_norm from the last iteration  */
    hapbnd = PetscMin(fgmres->haptol,hapbnd);
    if (tt > hapbnd) {
      /* scale new direction by its norm */
      ierr = VecScale(VEC_VV(loc_it+1),1.0/tt);CHKERRQ(ierr);
    } else {
      /* This happens when the solution is exactly reached. */
      /* So there is no new direction... */
      ierr   = VecSet(VEC_TEMP,0.0);CHKERRQ(ierr);     /* set VEC_TEMP to 0 */
      hapend = PETSC_TRUE;
    }
    /* note that for FGMRES we could get HES(loc_it+1, loc_it)  = 0 and the
       current solution would not be exact if HES was singular.  Note that
       HH non-singular implies that HES is no singular, and HES is guaranteed
       to be nonsingular when PREVECS are linearly independent and A is
       nonsingular (in GMRES, the nonsingularity of A implies the nonsingularity
       of HES). So we should really add a check to verify that HES is nonsingular.*/

    /* Now apply rotations to new col of hessenberg (and right side of system),
       calculate new rotation, and get new residual norm at the same time*/
    ierr = MyKSPFGMRESUpdateHessenberg(ksp,loc_it,hapend,&res_norm);CHKERRQ(ierr);
    if (ksp->reason) break;

    loc_it++;
    fgmres->it = (loc_it-1);   /* Add this here in case it has converged */

    ierr = PetscObjectSAWsTakeAccess((PetscObject)ksp);CHKERRQ(ierr);
    ksp->its++;
    ksp->rnorm = res_norm;
    ierr       = PetscObjectSAWsGrantAccess((PetscObject)ksp);CHKERRQ(ierr);

    ierr = (*ksp->converged)(ksp,ksp->its,res_norm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);

    /* Catch error in happy breakdown and signal convergence and break from loop */
    if (hapend) {
      if (!ksp->reason) {
        if (ksp->errorifnotconverged) SETERRQ1(PetscObjectComm((PetscObject)ksp),PETSC_ERR_NOT_CONVERGED,"You reached the happy break down, but convergence was not indicated. Residual norm = %g",(double)res_norm);
        else {
          ksp->reason = KSP_DIVERGED_BREAKDOWN;
          break;
        }
      }
    }
  }
  /* END OF ITERATION LOOP */
  ierr = KSPLogResidualHistory(ksp,res_norm);CHKERRQ(ierr);

  /*
     Monitor if we know that we will not return for a restart */
  if (loc_it && (ksp->reason || ksp->its >= ksp->max_it)) {
    ierr = KSPMonitor(ksp,ksp->its,res_norm);CHKERRQ(ierr);
  }

  if (itcount) *itcount = loc_it;

  /*
    Down here we have to solve for the "best" coefficients of the Krylov
    columns, add the solution values together, and possibly unwind the
    preconditioning from the solution
   */

  /* Form the solution (or the solution so far) */
  /* Note: must pass in (loc_it-1) for iteration count so that KSPFGMRESBuildSoln
     properly navigates */

  ierr = MyKSPFGMRESBuildSoln(RS(0),ksp->vec_sol,ksp->vec_sol,ksp,loc_it-1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESUpdateHessenberg"
PetscErrorCode MyKSPFGMRESUpdateHessenberg(KSP ksp,PetscInt it,PetscBool hapend,PetscReal *res)
{
  PetscScalar   *hh,*cc,*ss,tt;
  PetscInt      j;
  KSP_FGMRES    *fgmres = (KSP_FGMRES *)(ksp->data);

  PetscFunctionBegin;
  hh  = HH(0,it);  /* pointer to beginning of column to update - so
                      incrementing hh "steps down" the (it+1)th col of HH*/
  cc  = CC(0);     /* beginning of cosine rotations */
  ss  = SS(0);     /* beginning of sine rotations */

  /* Apply all the previously computed plane rotations to the new column
     of the Hessenberg matrix */
  /* Note: this uses the rotation [conj(c)  s ; -s   c], c= cos(theta), s= sin(theta),
     and some refs have [c   s ; -conj(s)  c] (don't be confused!) */

  for (j=1; j<=it; j++) {
    tt  = *hh;
#if defined(PETSC_USE_COMPLEX)
    *hh = PetscConj(*cc) * tt + *ss * *(hh+1);
#else
    *hh = *cc * tt + *ss * *(hh+1);
#endif
    hh++;
    *hh = *cc++ * *hh - (*ss++ * tt);
    /* hh, cc, and ss have all been incremented one by end of loop */
  }

  /*
    compute the new plane rotation, and apply it to:
     1) the right-hand-side of the Hessenberg system (RS)
        note: it affects RS(it) and RS(it+1)
     2) the new column of the Hessenberg matrix
        note: it affects HH(it,it) which is currently pointed to
        by hh and HH(it+1, it) (*(hh+1))
    thus obtaining the updated value of the residual...
  */

  /* compute new plane rotation */

  if (!hapend) {
#if defined(PETSC_USE_COMPLEX)
    tt        = PetscSqrtScalar(PetscConj(*hh) * *hh + PetscConj(*(hh+1)) * *(hh+1));
#else
    tt        = PetscSqrtScalar(*hh * *hh + *(hh+1) * *(hh+1));
#endif
    if (tt == 0.0) {
      ksp->reason = KSP_DIVERGED_NULL;
      PetscFunctionReturn(0);
    }

    *cc       = *hh / tt;   /* new cosine value */
    *ss       = *(hh+1) / tt;  /* new sine value */

    /* apply to 1) and 2) */
    *RS(it+1) = - (*ss * *RS(it));
#if defined(PETSC_USE_COMPLEX)
    *RS(it)   = PetscConj(*cc) * *RS(it);
    *hh       = PetscConj(*cc) * *hh + *ss * *(hh+1);
#else
    *RS(it)   = *cc * *RS(it);
    *hh       = *cc * *hh + *ss * *(hh+1);
#endif

    /* residual is the last element (it+1) of right-hand side! */
    *res      = PetscAbsScalar(*RS(it+1));

  } else { /* happy breakdown: HH(it+1, it) = 0, therfore we don't need to apply
            another rotation matrix (so RH doesn't change).  The new residual is
            always the new sine term times the residual from last time (RS(it)),
            but now the new sine rotation would be zero...so the residual should
            be zero...so we will multiply "zero" by the last residual.  This might
            not be exactly what we want to do here -could just return "zero". */

    *res = 0.0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESGetNewVectors"
PetscErrorCode MyKSPFGMRESGetNewVectors(KSP ksp,PetscInt it)
{
  KSP_FGMRES     *fgmres = (KSP_FGMRES *)ksp->data;
  PetscInt       nwork = fgmres->nwork_alloc; /* number of work vector chunks allocated */
  PetscInt       nalloc;                      /* number to allocate */
  PetscErrorCode ierr;
  PetscInt       k;

  PetscFunctionBegin;
  nalloc = fgmres->delta_allocate; /* number of vectors to allocate
                                      in a single chunk */
  /* Adjust the number to allocate to make sure that we don't exceed the
     number of available slots (fgmres->vecs_allocated)*/
  if (it + VEC_OFFSET + nalloc >= fgmres->vecs_allocated){
    nalloc = fgmres->vecs_allocated - it - VEC_OFFSET;
  }
  if (!nalloc) PetscFunctionReturn(0);

  fgmres->vv_allocated += nalloc; /* vv_allocated is the number of vectors allocated */

  /* work vectors */
  ierr = KSPCreateVecs(ksp,nalloc,&fgmres->user_work[nwork],0,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscLogObjectParents(ksp,nalloc,fgmres->user_work[nwork]);CHKERRQ(ierr);
  for (k=0; k < nalloc; k++) {
    fgmres->vecs[it+VEC_OFFSET+k] = fgmres->user_work[nwork][k];
  }
  /* specify size of chunk allocated */
  fgmres->mwork_alloc[nwork] = nalloc;

  /* preconditioned vectors */
  ierr = KSPCreateVecs(ksp,nalloc,&fgmres->prevecs_user_work[nwork],0,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscLogObjectParents(ksp,nalloc,fgmres->prevecs_user_work[nwork]);CHKERRQ(ierr);
  for (k=0; k < nalloc; k++) {
    fgmres->prevecs[it+VEC_OFFSET+k] = fgmres->prevecs_user_work[nwork][k];
  }

  /* increment the number of work vector chunks */
  fgmres->nwork_alloc++;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESBuildSoln"
PetscErrorCode MyKSPFGMRESBuildSoln(PetscScalar* nrs,Vec vguess,Vec vdest,KSP ksp,PetscInt it)
{
  PetscScalar    tt;
  PetscErrorCode ierr;
  PetscInt       ii,k,j;
  KSP_FGMRES     *fgmres = (KSP_FGMRES *)(ksp->data);

  PetscFunctionBegin;
  /* Solve for solution vector that minimizes the residual */

  /* If it is < 0, no fgmres steps have been performed */
  if (it < 0) {
    ierr = VecCopy(vguess,vdest);CHKERRQ(ierr); /* VecCopy() is smart, exists immediately if vguess == vdest */
    PetscFunctionReturn(0);
  }

  /* so fgmres steps HAVE been performed */

  /* solve the upper triangular system - RS is the right side and HH is
     the upper triangular matrix  - put soln in nrs */
  if (*HH(it,it) != 0.0) {
    nrs[it] = *RS(it) / *HH(it,it);
  } else {
    nrs[it] = 0.0;
  }
  for (ii=1; ii<=it; ii++) {
    k   = it - ii;
    tt  = *RS(k);
    for (j=k+1; j<=it; j++) tt  = tt - *HH(k,j) * nrs[j];
    nrs[k]   = tt / *HH(k,k);
  }

  /* Accumulate the correction to the soln of the preconditioned prob. in
     VEC_TEMP - note that we use the preconditioned vectors  */
  ierr = VecSet(VEC_TEMP,0.0);CHKERRQ(ierr); /* set VEC_TEMP components to 0 */
  ierr = VecMAXPY(VEC_TEMP,it+1,nrs,&PREVEC(0));CHKERRQ(ierr);

  /* put updated solution into vdest.*/
  if (vdest != vguess) {
    ierr = VecCopy(VEC_TEMP,vdest);CHKERRQ(ierr);
    ierr = VecAXPY(vdest,1.0,vguess);CHKERRQ(ierr);
  } else  {/* replace guess with solution */
    ierr = VecAXPY(vdest,1.0,VEC_TEMP);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
