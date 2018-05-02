#Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved#
ALL: all

<<<<<<< HEAD
#CFLAGS   =  -DEBUG  ${PETSC_CC_INCLUDES} -L/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/lib -I/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/include -lsmg2s2c
CFLAGS   =  ${PETSC_CC_INCLUDES} -L/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/lib -I/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/include -lsmg2s2c
SOURCEC  =  mat_generate.c mpi_lsa.c args_handler.c mpi_lsa_com.c read_matrix.c data_rw.c real2complex.c vec_generate.c 
SOURCEH  =  mat_generate.h mpi_lsa.h args_handler.h mpi_lsa_com.h read_matrix.h data_rw.h real2complex.h vec_generate.h 
CLEANFILE = mat_generate.o mpi_lsa.o args_handler.o mpi_lsa_com.o read_matrix.o data_rw.o real2complex.o vec_generate.o 
=======
CFLAGS   =  -DEBUG  ${PETSC_CC_INCLUDES}
SOURCEC  =  mpi_lsa.c args_handler.c mpi_lsa_com.c read_matrix.c data_rw.c real2complex.c vec_generate.c read_vector.c
SOURCEH  =  mpi_lsa.h args_handler.h mpi_lsa_com.h read_matrix.h data_rw.h real2complex.h vec_generate.h read_vector.h
CLEANFILE = mpi_lsa.o args_handler.o mpi_lsa_com.o read_matrix.o data_rw.o real2complex.o vec_generate.o read_vector.o
>>>>>>> e26a8fff0a5f33f8210759c863eb0da8b806ad37

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
