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

 #ifndef _READ_MATRIX_H
#define _READ_MATRIX_H

#include "petscmat.h"
#include "petscerror.h"

#endif
PetscErrorCode generate_canoic_vector(PetscInt size, PetscScalar value, PetscInt position, Vec * v);
PetscErrorCode generate_homogene_vector(PetscInt size, PetscInt scal,Vec * v);
PetscErrorCode generate_random_vector(PetscInt size, PetscInt low, PetscInt up,Vec * v);
PetscErrorCode generate_continuzero_vector(PetscInt size, PetscInt start, PetscInt length, PetscScalar value, Vec * v);
PetscErrorCode generate_random_seed_vector(PetscInt size, PetscReal low, PetscReal up, PetscReal seed, Vec * v);
