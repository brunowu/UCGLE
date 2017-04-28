
#This file is part of software for the implementation of UCGLE method, under the supervision of Serge G. Petiton
#<serge.petiton@univ-lille1.fr>.
 
#Copyright (C) 2011—. Pierre-Yves AQUILANTI and Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr> in Maison de la Simulation. 
#All rights reserved.
 
#Permission to use, copy, modify and distribute this software for personal and educational use is hereby granted
#without fee, provided that the above copyright notice appears in all copies and that both that copyright notice 
#and this permission notice appear in supporting documentation, and that the names of all authors are not used in 
#advertising or publicity pertaining to distribution of the software without specific, written prior permission. 
#Xinzhe WU and the author make no representations about the suitability of this software for any purpose. It is 
#provided "as is" without express or implied warranty.
 
#You should have received a copy of the GNU Lesser General Public License along with UCGLE.  If not, see 
#<http://www.gnu.org/licenses/>.

#For more information, contact with Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr>.
 
 

 ALL: all

CFLAGS   =  -DEBUG  ${PETSC_CC_INCLUDES}
SOURCEC  =  mpi_lsa.c args_handler.c mpi_lsa_com.c read_matrix.c data_rw.c real2complex.c vec_generate.c
SOURCEH  =  mpi_lsa.h args_handler.h mpi_lsa_com.h read_matrix.h data_rw.h real2complex.h vec_generate.h
CLEANFILE = mpi_lsa.o args_handler.o mpi_lsa_com.o read_matrix.o data_rw.o real2complex.o vec_generate.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
