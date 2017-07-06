#Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved#
ALL: all

CFLAGS   = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -DDEBUG
SOURCEC  =  gmres_cycle.c gmres_solve.c kspsolve.c gmres.c gmres_precond.c classical_gmres.c
SOURCEH  =  gmres_cycle.h gmres_solve.h kspsolve.h gmres.h gmres_precond.h classical_gmres.h
CLEANFILE    = gmres_cycle.o gmres_solve.o kspsolve.o gmres.o gmres_precond.o classical_gmres.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
