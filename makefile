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
CFLAGS	= -O0
#-I${PETSC_DIR}/src/ksp/ksp/impls/gmres/fgmres/

CLEANFILES  = ${EXEC}
HFILES = -I./Libs  -I./Solvers/Utils -I./Solvers/Arnoldi -I./Tuning
OFILES= ${wildcard ./Libs/*.o} ${wildcard ./Solvers/*/*.o} ${wildcard ./Tuning/*.o}

#matrices directory
MDIR=./data

##################################################################
##################      Solvers FLAGS      #######################
##################################################################

#debug options

#DEBUG_VALGRIND = valgrind --tool=memcheck -q
#DEBUG_KSP_VIEW = -ksp_view

RESTART_MAX = 20
GMRES_PRECISION= 1e-10
GMRES_RESTART= ${RESTART_MAX}
GMRES_NB_NODES=1
GMRES_MONITOR= -ksp_monitor_true_residual
NTIMES = 1
GMRES_FLAGS= -ksp_rtol 1e-100 -ksp_divtol 1e1000 -ksp_max_it 20000 -pc_type none -ksp_atol ${GMRES_PRECISION} -ksp_gmres_restart ${GMRES_RESTART}\
		${GMRES_MONITOR} -lsa_gmres ${GMRES_NB_NODES} -ntimes ${NTIMES} \
		#-log_summary

#arnoldi options
ARNOLDI_PRECISION= 1e-1
ARNOLDI_NBEIGEN= 18
ARNOLDI_NB_NODES=1
ARNOLDI_MONITOR = -eps_monitor_conv
ARNOLDI_FLAGS= -eps_ncv 20 -eps_type arnoldi -eps_true_residual -eps_largest_imaginary -eps_nev ${ARNOLDI_NBEIGEN} -eps_tol ${ARNOLDI_PRECISION} \
                ${ARNOLDI_MONITOR} -lsa_arnoldi ${ARNOLDI_NB_NODES} -eps_max_it 50
#ls options
LS_POWER = 10
LS_POLY_APPL = 5
LS_LATENCY=1
LS_PC_USE =0
LS_NO_USE_LS= -ksp_ls_nols
LS_HANG_IT= 20000
LS_HANG_TIME=  100000
# LS_LOAD_ANY = -ksp_ls_load_any

#LS_FLAGS = -ksp_ls_power ${LS_POWER} ${LS_NO_USE_LS}-ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY} -ksp_ls_cexport

#LS_FLAGS = -ksp_ls_power ${LS_POWER} ${LS_NO_USE_LS}-ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY} -ksp_ls_load ./lsqr.bin

LS_FLAGS = -ksp_ls_power ${LS_POWER} ${LS_NO_USE_LS}-ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY}

#final flag composition  ${DEBUGG}

GLSA_FLAGS = ${GMRES_FLAGS} ${ARNOLDI_FLAGS} ${LS_FLAGS} ${DEBUG_KSP_VIEW} ${RUN_FLAGS}
MPI_NODES = ${shell echo ${GMRES_NB_NODES}+${ARNOLDI_NB_NODES}+2 | bc}



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

runutm300:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/utm300.mtx_300x300_3155nnz

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
        -mfile ${MDIR}/Eigen_known_matrix_nb_300_300x300_10264_nnz
