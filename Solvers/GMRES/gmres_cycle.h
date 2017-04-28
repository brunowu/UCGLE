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

#ifndef GMRES_CYCLE_H
#define GMRES_CYCLE_H

/* just one include necessary */
#include "fgmresimpl.h"
#include "petsc.h"
#include "petsc/private/kspimpl.h"


PetscErrorCode MyKSPFGMRESCycle(PetscInt *itcount,KSP ksp);
PetscErrorCode MyKSPFGMRESGetNewVectors(KSP,PetscInt);
PetscErrorCode MyKSPFGMRESUpdateHessenberg(KSP,PetscInt,PetscBool,PetscReal *);
PetscErrorCode MyKSPFGMRESBuildSoln(PetscScalar*,Vec,Vec,KSP,PetscInt);

#endif
