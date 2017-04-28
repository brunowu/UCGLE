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

#include "gmres_solve.h"

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESResidual"
static PetscErrorCode MyKSPFGMRESResidual(KSP ksp)
{
  KSP_FGMRES     *fgmres = (KSP_FGMRES*)(ksp->data);
  Mat            Amat,Pmat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);

  /* put A*x into VEC_TEMP */
  ierr = KSP_MatMult(ksp,Amat,ksp->vec_sol,VEC_TEMP);CHKERRQ(ierr);
  /* now put residual (-A*x + f) into vec_vv(0) */
  ierr = VecWAXPY(VEC_VV(0),-1.0,VEC_TEMP,ksp->vec_rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyKSPSolve_FGMRES"
PetscErrorCode MyKSPSolve_FGMRES(KSP ksp,com_lsa * com)
{
  PetscErrorCode ierr;
  PetscInt       cycle_its = 0; /* iterations done in a call to KSPFGMRESCycle */
  KSP_FGMRES     *fgmres   = (KSP_FGMRES*)ksp->data;
  PetscBool      diagonalscale;
  Vec vec_tmp;

  PetscFunctionBegin;
  ierr = PCGetDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
  if (diagonalscale) SETERRQ1(PetscObjectComm((PetscObject)ksp),PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);

  ierr     = PetscObjectSAWsTakeAccess((PetscObject)ksp);CHKERRQ(ierr);
  ksp->its = 0;
  ierr     = PetscObjectSAWsGrantAccess((PetscObject)ksp);CHKERRQ(ierr);

  /* Compute the initial (NOT preconditioned) residual */
  if (!ksp->guess_zero) {
    ierr = MyKSPFGMRESResidual(ksp);CHKERRQ(ierr);
  } else { /* guess is 0 so residual is F (which is in ksp->vec_rhs) */
    ierr = VecCopy(ksp->vec_rhs,VEC_VV(0));CHKERRQ(ierr);
  }
  /* now the residual is in VEC_VV(0) - which is what
     KSPFGMRESCycle expects... */

  ierr = MyKSPFGMRESCycle(&cycle_its,ksp);CHKERRQ(ierr);
  while (!ksp->reason) {
    if(!GmresLSAPrecond(com,ksp))
    {
///	    PetscPrintf(PETSC_COMM_WORLD,"--->>>Preconditioning of LS method in: %d iterations\n\n",ksp->its);
    }

    ierr = MyKSPFGMRESResidual(ksp);CHKERRQ(ierr);
    if (ksp->its >= ksp->max_it) break;
    ierr = MyKSPFGMRESCycle(&cycle_its,ksp);CHKERRQ(ierr);
    ierr = VecDuplicate(ksp->vec_rhs,&vec_tmp);CHKERRQ(ierr);
    mpi_lsa_com_vec_send(com,&vec_tmp);
  }

  /* mark lack of convergence */
  if (ksp->its >= ksp->max_it && !ksp->reason) ksp->reason = KSP_DIVERGED_ITS;
  PetscFunctionReturn(0);
}
