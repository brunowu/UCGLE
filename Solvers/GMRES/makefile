#Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved#
ALL: all gmres
GMRES_EXEC: GMRES

CFLAGS   = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -DDEBUG
SOURCEC  =  classical_gmres.c
SOURCEH  =  classical_gmres.h
CLEANFILE    = classical_gmres.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

gmres: classical_gmres.o
	-@echo "BEGINNING TO COMPILE SLAVE1 "
	-@echo "========================================="
	-@${CLINKER} -o ${GMRES_EXEC} classical_gmres.o ${PETSC_LIB}
	-@echo "Completed SLAVE1 Compilation"
	-@echo "========================================="