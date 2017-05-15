#Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved
ALL: all

CFLAGS   = -DEBUG -I${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/${PETSC_ARCH}/lib/
SOURCEC  =  father.c
SOURCEH  =  father.h
CLEANFILE  = father.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
