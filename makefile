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
##################      Solvers FLAGS      #######################
##################################################################

#################       Debug flags    ########################

#DEBUG_VALGRIND = valgrind --tool=memcheck -q
#DEBUG_KSP_VIEW = -ksp_view

#################       General Flags    ########################
#MAT = EBMG_matrix_nb_200000_200000x200000_1.41977e+07_nnz
#MAT = EBMG_matrix_nb_400000_400000x400000_2.83976e+07_nnz
#MAT = EBMG_matrix_nb_800000_800000x800000_5.67979e+07_nnz
#MAT = EBMG_matrix_nb_1600000_1600000x1600000_1.13598e+08_nnz
#MAT = Eigen_known_matrix_nb_500_500x500_40425_nnz
#MAT = EDMG_matrix_nb_500_500x500_40419_nnz
#MAT = EBMG_matrix_nb_500_500x500_33039_nnz
#MAT = EBMG_matrix_nb_500_500x500_33063_nnz
#MAT = EBMG_matrix_nb_500_500x500_33098_nnz
#MAT = EBMG_matrix_nb_500_500x500_33036_nnz
#MAT = matline_nb_6_1800x1800_19216_nnz
#MAT = matline_nb_6_1800x1800_18933_nnzK
#MAT = utm300_300x300.dat
MAT = utm300.mtx_300x300_3155nnz
#MAT = matline_nb_3_900x900_9468_nnz
#MAT = matblock_nb_3_900x900_15775_nnz

NTIMES = 1
INITIAL_GUESS_SEQ_RHS = -initial_guess_seq_rhs
INITIAL_GUESS_POWER = 5


#################       MPI Flags        ########################

GMRES_NB_NODES=1
ARNOLDI_NB_NODES=1

MPI_NODES = ${shell echo ${GMRES_NB_NODES}+${ARNOLDI_NB_NODES}+ 2 | bc}
LSA_GMRES =-lsa_gmres
LSA_ARNOLDI=-lsa_arnoldi

#classical GMRES MPI flags
#MPI_NODES = 2
#################       GMRES Flags       ########################

RESTART_MAX =150
GMRES_PRECISION= 1e-10
GMRES_RESTART= ${RESTART_MAX}
GMRES_MONITOR= -ksp_monitor_true_residual
KSP_MAX_ITS = 20000
PC_TYPE = none
#PC_TYPE=jacobi -pc_jacobi_type rowmax -pc_jacobi_abs
#PC_TYPE = asm -pc_asm_block 2 -pc_asm_type none
#PC_TYPE = sor -pc_sor_its 5 -pc_sor_omega 0.179
#CUDA_TYPE = -mat_type aijcusparse -vec_type cuda
#GMRES_FT = -GMRES_FT
GMRES_FLAGS= -ksp_rtol 1e-100 -ksp_divtol 1e1000 -ksp_max_it ${KSP_MAX_ITS} -pc_type ${PC_TYPE} -ksp_atol ${GMRES_PRECISION} -ksp_gmres_restart ${GMRES_RESTART}\
		${GMRES_MONITOR} ${LSA_GMRES} ${GMRES_NB_NODES} -ntimes ${NTIMES} ${CUDA_TYPE}\
		${GMRES_FT} ${INITIAL_GUESS_SEQ_RHS} -initial_guess_power ${INITIAL_GUESS_POWER}

#		-log_view

#################       ERAM Flags         ########################

ARNOLDI_PRECISION= 1e-1
ARNOLDI_NBEIGEN=3
#ARNOLDI_MONITOR = -eps_monitor_conv
#ARNOLDI_FT_SIM = -ArnoldiFT 4

ARNOLDI_FLAGS= -eps_ncv 60 -eps_type arnoldi -eps_true_residual -eps_largest_imaginary -eps_nev ${ARNOLDI_NBEIGEN} -eps_tol ${ARNOLDI_PRECISION} \
                ${ARNOLDI_MONITOR} ${LSA_ARNOLDI} ${ARNOLDI_NB_NODES} -eps_max_it 50 \
		${ARNOLDI_FT_SIM} 
#################       LS Flags           ########################

LS_POWER = 10
LS_POLY_APPL = 10
LS_LATENCY=1
LS_PC_USE =1
LS_NO_USE_LS= -ksp_ls_nols
LS_HANG_IT= 20000
LS_HANG_TIME=  100000
#################    LS Reusability Flags  ########################

#LS_CEXPORT = -ksp_ls_cexport
#LS_LOAD = -ksp_ls_load
#LS_LOAD_FILE = ./lsqr.bin

LS_FLAGS = -ksp_ls_power ${LS_POWER} ${LS_NO_USE_LS}-ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY} ${LS_CEXPORT} ${LS_LOAD} ${LS_LOAD_FILE}
#################      Flags Combinations      #####################

GLSA_FLAGS = ${GMRES_FLAGS} ${ARNOLDI_FLAGS} ${LS_FLAGS} ${DEBUG_KSP_VIEW} ${RUN_FLAGS}

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
	#@${CLINKER} -I${SLEPC_DIR}/include ${OFILES} ${HFILES} -L${SLEPC_LIB} -L. -L${SLEPC_DIR}/${PETSC_ARCH}/lib -g -v  -o convert convertor.o
	#@${CLINKER} -L{PETSC_DIR}/${PETSC_ARCH}/lib -g -v -o ${EXEC} main.o ${OFILES} ${HFILES} -I${PETSC_DIR}/include -lm
	#${CLINKER}  -L${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/ -L. ${OFILES}  -o ${EXEC} main.o ${HFILES} -I${PETSC_DIR}/include/petsc/private/ -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include

	-@echo "Completed building application"
	-@echo "========================================="



effacer :
	-rm *.bin
	-rm *.o
	-rm ./*/*.o
	-rm ./*/*/*.o


	-@echo "Completed building application"
	-@echo "========================================="

##################################################################
##################     Execution Rules     #######################
##################################################################

runexec:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/${MAT}

runmatBlock:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ${MDIR}/matBlock_nb_300_90000x90000_1.88984e+06_nnz

runmatLine:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ${MDIR}/matLine_nb_600_180000x180000_1.91215e+06_nnz \
        2>&1 | tee log.txt

runa:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ${MDIR}/Eigen_known_matrix_nb_274_274x274_3179_nnz

runb:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ${MDIR}/Eigen_known_matrix_nb_274_274x274_5734_nnz
runc:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ${MDIR}/EBMG_matrix_nb_100_100x100_1498_nnz
