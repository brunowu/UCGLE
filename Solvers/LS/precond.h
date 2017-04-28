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

#include "petsc.h"
#include "petscpc.h"
#include "petscmat.h"
#include "petscksp.h"
#include <stdlib.h>
#include "../Utils/lib.h"

#define at(a, i, j) (a)[j][i]

PetscErrorCode LSPrecond(PetscReal a_ell, PetscReal d_ell, PetscReal c_ell,
			PetscScalar * eta, PetscScalar * alpha, PetscScalar * beta,
			PetscScalar * delta, PetscScalar * c, PetscScalar * d,
			PetscInt * mu, PetscInt * nb_eigen, PetscInt * min_eigen, PetscInt * nb_eigen_all);
