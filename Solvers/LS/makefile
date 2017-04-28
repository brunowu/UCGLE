#Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved
ALL: all

CFLAGS   =  -DDEBUG -I${PETSC_DIR}/include ${PETSC_LIB} -I${PETSC_CC_INCLUDES}
SOURCEC  =  lsqr.c precond.c
SOURCEH  =  lsqr.h precond.h
CLEANFILE    = lsqr.o precond.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
