#Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved#
ALL: all

CFLAGS   = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -DDEBUG
SOURCEC  =  arnoldi.c
SOURCEH  =  arnoldi.h
CLEANFILE    = arnoldi.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
