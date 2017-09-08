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

ALL:  blib exec

#compilation flags
DIRS    = Libs Tuning Solvers
EXEC    = hyperh
CFLAGS	= -O0 -O3
#-I${PETSC_DIR}/src/ksp/ksp/impls/gmres/fgmres/

CLEANFILES  = ${EXEC}
HFILES = -I./Libs  -I./Solvers/Utils -I./Solvers/Arnoldi -I./Tuning
OFILES= ${wildcard ./Libs/*.o} ${wildcard ./Solvers/*/*.o} ${wildcard ./Tuning/*.o}

#matrices directory
MDIR=./data

##################################################################
##################   Compilation rules     #######################
##################################################################

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

blib:
	@for a in $(DIRS); do \
        if [ -d $$a ]; then \
        echo "processing folder $$a"; \
        $(MAKE) -C $$a; \
        fi; \
        done;
	@echo "$$a Done!"
	-@echo "Completed building libraries"
	-@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES} main.o
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="

rmat :
	-rm -drf matrices_tmp/*


exec: main.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	@${CLINKER} -I${SLEPC_DIR}/include ${OFILES} ${HFILES} -L${SLEPC_LIB} -L. -L${SLEPC_DIR}/${PETSC_ARCH}/lib -g -v -o ${EXEC} main.o
	-@echo "Completed building application"
	-@echo "========================================="

effacer :
	-rm *.bin
	-rm *.o
	-rm ./*/*.o
	-rm ./*/*/*.o
	-@echo "Completed building application"
	-@echo "========================================="
