#Copyright (c) 2011—2016. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved
ALL: all

CFLAGS   = ${PETSC_CC_INCLUDES}
SOURCEC  =  ellipse.c convhul.c triep.c lib.c
SOURCEH  =  ellipse.h convhul.h triep.h lib.h
CLEANFILE    = ellipse.o convhul.o triep.o lib.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@
all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
