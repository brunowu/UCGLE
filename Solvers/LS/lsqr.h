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

#ifndef LSQR_H
#define LSQR_H 

#include "petsc.h"
#include "../../Libs/mpi_lsa.h"
#include "../../Libs/mpi_lsa_com.h"
#include "../../Libs/data_rw.h"
#include "../Utils/triep.h"
#include "../Utils/lib.h"
#include "../Utils/convhul.h"
#include "../Utils/ellipse.h"
#include "./precond.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef EIGEN_MIN
#define EIGEN_MIN 5
#endif

#ifndef EIGEN_ALL
#define EIGEN_ALL 20	
#endif

#ifndef EIGEN_MAX
#define EIGEN_MAX 20	
#endif

PetscErrorCode LSQR(com_lsa * com, int * vector_size);

#endif
