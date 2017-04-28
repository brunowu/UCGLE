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

#include "kspsolve.h"
#undef __FUNCT__
#define __FUNCT__ "MyKSPSolve"
PetscErrorCode MyKSPSolve(KSP ksp,Vec b,Vec x,com_lsa * com)
{
  PetscErrorCode    ierr;
  PetscMPIInt       rank;
  PetscBool         flag1,flag2,flag3,flg = PETSC_FALSE,inXisinB=PETSC_FALSE,guess_zero;
  Mat               mat,premat;
  MPI_Comm          comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  if (b) PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  if (x) PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  comm = PetscObjectComm((PetscObject)ksp);
  if (x && x == b) {
    if (!ksp->guess_zero) SETERRQ(comm,PETSC_ERR_ARG_INCOMP,"Cannot use x == b with nonzero initial guess");
    ierr     = VecDuplicate(b,&x);CHKERRQ(ierr);
    inXisinB = PETSC_TRUE;
  }
  if (b) {
    ierr         = PetscObjectReference((PetscObject)b);CHKERRQ(ierr);
    ierr         = VecDestroy(&ksp->vec_rhs);CHKERRQ(ierr);
    ksp->vec_rhs = b;
  }
  if (x) {
    ierr         = PetscObjectReference((PetscObject)x);CHKERRQ(ierr);
    ierr         = VecDestroy(&ksp->vec_sol);CHKERRQ(ierr);
    ksp->vec_sol = x;
  }
  ierr = KSPViewFromOptions(ksp,NULL,"-ksp_view_pre");CHKERRQ(ierr);

  if (ksp->presolve) {
    ierr = (*ksp->presolve)(ksp,ksp->vec_rhs,ksp->vec_sol,ksp->prectx);CHKERRQ(ierr);
  }
  ierr = PetscLogEventBegin(KSP_Solve,ksp,ksp->vec_rhs,ksp->vec_sol,0);CHKERRQ(ierr);

  /* reset the residual history list if requested */
  if (ksp->res_hist_reset) ksp->res_hist_len = 0;
  ksp->transpose_solve = PETSC_FALSE;

  if (ksp->guess) {
    ierr            = KSPFischerGuessFormGuess(ksp->guess,ksp->vec_rhs,ksp->vec_sol);CHKERRQ(ierr);
    ksp->guess_zero = PETSC_FALSE;
  }
  /* KSPSetUp() scales the matrix if needed */
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPSetUpOnBlocks(ksp);CHKERRQ(ierr);
  VecLocked(ksp->vec_sol,3);

  /* diagonal scale RHS if called for */
  if (ksp->dscale) {

    ierr = VecPointwiseMult(ksp->vec_rhs,ksp->vec_rhs,ksp->diagonal);CHKERRQ(ierr);
    /* second time in, but matrix was scaled back to original */

    if (ksp->dscalefix && ksp->dscalefix2) {
      Mat mat,pmat;

      ierr = PCGetOperators(ksp->pc,&mat,&pmat);CHKERRQ(ierr);
      ierr = MatDiagonalScale(pmat,ksp->diagonal,ksp->diagonal);CHKERRQ(ierr);
      if (mat != pmat) {ierr = MatDiagonalScale(mat,ksp->diagonal,ksp->diagonal);CHKERRQ(ierr);}
    }

    /*  scale initial guess */
    if (!ksp->guess_zero) {
      if (!ksp->truediagonal) {
        ierr = VecDuplicate(ksp->diagonal,&ksp->truediagonal);CHKERRQ(ierr);
        ierr = VecCopy(ksp->diagonal,ksp->truediagonal);CHKERRQ(ierr);
        ierr = VecReciprocal(ksp->truediagonal);CHKERRQ(ierr);
      }
      ierr = VecPointwiseMult(ksp->vec_sol,ksp->vec_sol,ksp->truediagonal);CHKERRQ(ierr);
    }
  }

  ierr = PCPreSolve(ksp->pc,ksp);CHKERRQ(ierr);

  if (ksp->guess_zero) { ierr = VecSet(ksp->vec_sol,0.0);CHKERRQ(ierr);}
  if (ksp->guess_knoll) {
    ierr            = PCApply(ksp->pc,ksp->vec_rhs,ksp->vec_sol);CHKERRQ(ierr);
    ierr            = KSP_RemoveNullSpace(ksp,ksp->vec_sol);CHKERRQ(ierr);
    ksp->guess_zero = PETSC_FALSE;
  }

  /* can we mark the initial guess as zero for this solve? */
  guess_zero = ksp->guess_zero;
  if (!ksp->guess_zero) {
    PetscReal norm;

    ierr = VecNormAvailable(ksp->vec_sol,NORM_2,&flg,&norm);CHKERRQ(ierr);
    if (flg && !norm) ksp->guess_zero = PETSC_TRUE;
  }

  ierr = VecLockPush(ksp->vec_rhs);CHKERRQ(ierr);

/* !!!!!!!!!!!!!!!!!!!!!!    IMPORTANT   !!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!! THIS IS A MAJOR CHANGE OF THE ORIGINAL SOURCE CODE AND IT MAY ENGENDER SOME PROBLEMS  <<<< So keep an eye on it >>>> !!!!!!!!!
this is the default way to launch the actual solver:
//  ierr = (*ksp->ops->solve)(ksp);CHKERRQ(ierr);
 the problem is that we set KSPSetType to fgmres what makes PETSC use its own solver
So in order to use our FGMRES we make it explicitly like follows
*/

  ierr = MyKSPSolve_FGMRES(ksp,com);CHKERRQ(ierr);
  ierr = VecLockPop(ksp->vec_rhs);CHKERRQ(ierr);
  ksp->guess_zero = guess_zero;

  if (!ksp->reason) SETERRQ(comm,PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");
  ierr = KSPReasonViewFromOptions(ksp);CHKERRQ(ierr);
  ierr = PCPostSolve(ksp->pc,ksp);CHKERRQ(ierr);

  /* diagonal scale solution if called for */
  if (ksp->dscale) {
    ierr = VecPointwiseMult(ksp->vec_sol,ksp->vec_sol,ksp->diagonal);CHKERRQ(ierr);
    /* unscale right hand side and matrix */
    if (ksp->dscalefix) {
      Mat mat,pmat;

      ierr = VecReciprocal(ksp->diagonal);CHKERRQ(ierr);
      ierr = VecPointwiseMult(ksp->vec_rhs,ksp->vec_rhs,ksp->diagonal);CHKERRQ(ierr);
      ierr = PCGetOperators(ksp->pc,&mat,&pmat);CHKERRQ(ierr);
      ierr = MatDiagonalScale(pmat,ksp->diagonal,ksp->diagonal);CHKERRQ(ierr);
      if (mat != pmat) {ierr = MatDiagonalScale(mat,ksp->diagonal,ksp->diagonal);CHKERRQ(ierr);}
      ierr            = VecReciprocal(ksp->diagonal);CHKERRQ(ierr);
      ksp->dscalefix2 = PETSC_TRUE;
    }
  }
  ierr = PetscLogEventEnd(KSP_Solve,ksp,ksp->vec_rhs,ksp->vec_sol,0);CHKERRQ(ierr);
  if (ksp->postsolve) {
    ierr = (*ksp->postsolve)(ksp,ksp->vec_rhs,ksp->vec_sol,ksp->postctx);CHKERRQ(ierr);
  }

  if (ksp->guess) {
    ierr = KSPFischerGuessUpdate(ksp->guess,ksp->vec_sol);CHKERRQ(ierr);
  }

  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  ierr = PCGetOperators(ksp->pc,&mat,&premat);CHKERRQ(ierr);
  ierr = MatViewFromOptions(mat,(PetscObject)ksp,"-ksp_view_mat");CHKERRQ(ierr);
  ierr = MatViewFromOptions(premat,(PetscObject)ksp,"-ksp_view_pmat");CHKERRQ(ierr);
  ierr = VecViewFromOptions(ksp->vec_rhs,(PetscObject)ksp,"-ksp_view_rhs");CHKERRQ(ierr);

  flag1 = PETSC_FALSE;
  flag2 = PETSC_FALSE;
  flag3 = PETSC_FALSE;
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_compute_eigenvalues",&flag1,NULL);CHKERRQ(ierr);
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_plot_eigenvalues",&flag2,NULL);CHKERRQ(ierr);
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_plot_eigencontours",&flag3,NULL);CHKERRQ(ierr);
  if (flag1 || flag2 || flag3) {
    PetscInt  nits,n,i,neig;
    PetscReal *r,*c;

    ierr = KSPGetIterationNumber(ksp,&nits);CHKERRQ(ierr);
    n    = nits+2;

    if (!nits) {
      ierr = PetscPrintf(comm,"Zero iterations in solver, cannot approximate any eigenvalues\n");CHKERRQ(ierr);
    } else {
      ierr = PetscMalloc2(n,&r,n,&c);CHKERRQ(ierr);
      ierr = KSPComputeEigenvalues(ksp,n,r,c,&neig);CHKERRQ(ierr);
      if (flag1) {
        ierr = PetscPrintf(comm,"Iteratively computed eigenvalues\n");CHKERRQ(ierr);
        for (i=0; i<neig; i++) {
          if (c[i] >= 0.0) {
            ierr = PetscPrintf(comm,"%g + %gi\n",(double)r[i],(double)c[i]);CHKERRQ(ierr);
          } else {
            ierr = PetscPrintf(comm,"%g - %gi\n",(double)r[i],-(double)c[i]);CHKERRQ(ierr);
          }
        }
      }
      if (flag2 && !rank) {
        PetscDraw   draw;
        PetscDrawSP drawsp;

        if (!ksp->eigviewer) {
          ierr = PetscViewerDrawOpen(PETSC_COMM_SELF,0,"Iteratively Computed Eigenvalues",PETSC_DECIDE,PETSC_DECIDE,400,400,&ksp->eigviewer);CHKERRQ(ierr);
        }
        ierr = PetscViewerDrawGetDraw(ksp->eigviewer,0,&draw);CHKERRQ(ierr);
        ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
        ierr = PetscDrawSPReset(drawsp);CHKERRQ(ierr);
        for (i=0; i<neig; i++) {
          ierr = PetscDrawSPAddPoint(drawsp,r+i,c+i);CHKERRQ(ierr);
        }
        ierr = PetscDrawSPDraw(drawsp,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
      }
      if (flag3 && !rank) {
        ierr = KSPPlotEigenContours_Private(ksp,neig,r,c);CHKERRQ(ierr);
      }
      ierr = PetscFree2(r,c);CHKERRQ(ierr);
    }
  }

  flag1 = PETSC_FALSE;
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_compute_singularvalues",&flag1,NULL);CHKERRQ(ierr);
  if (flag1) {
    PetscInt nits;

    ierr = KSPGetIterationNumber(ksp,&nits);CHKERRQ(ierr);
    if (!nits) {
      ierr = PetscPrintf(comm,"Zero iterations in solver, cannot approximate any singular values\n");CHKERRQ(ierr);
    } else {
      PetscReal emax,emin;

      ierr = KSPComputeExtremeSingularValues(ksp,&emax,&emin);CHKERRQ(ierr);
      ierr = PetscPrintf(comm,"Iteratively computed extreme singular values: max %g min %g max/min %g\n",(double)emax,(double)emin,(double)(emax/emin));CHKERRQ(ierr);
    }
  }


  flag1 = PETSC_FALSE;
  flag2 = PETSC_FALSE;
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_compute_eigenvalues_explicitly",&flag1,NULL);CHKERRQ(ierr);
  ierr  = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_plot_eigenvalues_explicitly",&flag2,NULL);CHKERRQ(ierr);
  if (flag1 || flag2) {
    PetscInt  n,i;
    PetscReal *r,*c;
    ierr = VecGetSize(ksp->vec_sol,&n);CHKERRQ(ierr);
    ierr = PetscMalloc2(n,&r,n,&c);CHKERRQ(ierr);
    ierr = KSPComputeEigenvaluesExplicitly(ksp,n,r,c);CHKERRQ(ierr);
    if (flag1) {
      ierr = PetscPrintf(comm,"Explicitly computed eigenvalues\n");CHKERRQ(ierr);
      for (i=0; i<n; i++) {
        if (c[i] >= 0.0) {
          ierr = PetscPrintf(comm,"%g + %gi\n",(double)r[i],(double)c[i]);CHKERRQ(ierr);
        } else {
          ierr = PetscPrintf(comm,"%g - %gi\n",(double)r[i],-(double)c[i]);CHKERRQ(ierr);
        }
      }
    }
    if (flag2 && !rank) {
      PetscDraw   draw;
      PetscDrawSP drawsp;

      if (!ksp->eigviewer) {
        ierr = PetscViewerDrawOpen(PETSC_COMM_SELF,0,"Explicitly Computed Eigenvalues",0,320,400,400,&ksp->eigviewer);CHKERRQ(ierr);
      }
      ierr = PetscViewerDrawGetDraw(ksp->eigviewer,0,&draw);CHKERRQ(ierr);
      ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
      ierr = PetscDrawSPReset(drawsp);CHKERRQ(ierr);
      for (i=0; i<n; i++) {
        ierr = PetscDrawSPAddPoint(drawsp,r+i,c+i);CHKERRQ(ierr);
      }
      ierr = PetscDrawSPDraw(drawsp,PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
    }
    ierr = PetscFree2(r,c);CHKERRQ(ierr);
  }

  ierr = PetscOptionsHasName(((PetscObject)ksp)->prefix,"-ksp_view_mat_explicit",&flag2);CHKERRQ(ierr);
  if (flag2) {
    Mat A,B;
    ierr = PCGetOperators(ksp->pc,&A,NULL);CHKERRQ(ierr);
    ierr = MatComputeExplicitOperator(A,&B);CHKERRQ(ierr);
    ierr = MatViewFromOptions(B,(PetscObject)ksp,"-ksp_view_mat_explicit");CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
  }
  ierr = PetscOptionsHasName(((PetscObject)ksp)->prefix,"-ksp_view_preconditioned_operator_explicit",&flag2);CHKERRQ(ierr);
  if (flag2) {
    Mat B;
    ierr = KSPComputeExplicitOperator(ksp,&B);CHKERRQ(ierr);
    ierr = MatViewFromOptions(B,(PetscObject)ksp,"-ksp_view_preconditioned_operator_explicit");CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
  }
  ierr = KSPViewFromOptions(ksp,NULL,"-ksp_view");CHKERRQ(ierr);

  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(((PetscObject)ksp)->prefix,"-ksp_final_residual",&flg,NULL);CHKERRQ(ierr);
  if (flg) {
    Mat       A;
    Vec       t;
    PetscReal norm;
    if (ksp->dscale && !ksp->dscalefix) SETERRQ(comm,PETSC_ERR_ARG_WRONGSTATE,"Cannot compute final scale with -ksp_diagonal_scale except also with -ksp_diagonal_scale_fix");
    ierr = PCGetOperators(ksp->pc,&A,NULL);CHKERRQ(ierr);
    ierr = VecDuplicate(ksp->vec_rhs,&t);CHKERRQ(ierr);
    ierr = KSP_MatMult(ksp,A,ksp->vec_sol,t);CHKERRQ(ierr);
    ierr = VecAYPX(t, -1.0, ksp->vec_rhs);CHKERRQ(ierr);
    ierr = VecNorm(t,NORM_2,&norm);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
    ierr = PetscPrintf(comm,"KSP final norm of residual %g\n",(double)norm);CHKERRQ(ierr);
  }
  ierr = VecViewFromOptions(ksp->vec_sol,(PetscObject)ksp,"-ksp_view_solution");CHKERRQ(ierr);

  if (inXisinB) {
    ierr = VecCopy(x,b);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
  }
  ierr = PetscObjectSAWsBlock((PetscObject)ksp);CHKERRQ(ierr);
  if (ksp->errorifnotconverged && ksp->reason < 0) SETERRQ(comm,PETSC_ERR_NOT_CONVERGED,"KSPSolve has not converged");
  PetscFunctionReturn(0);
}
