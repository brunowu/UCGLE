#!/bin/bash

EXEC=hyperh

##################################################################
##################      Matrix Location      #####################
##################################################################

MDIR=./data

##################################################################
##################      Solvers FLAGS      #######################
##################################################################

MAT = utm300.mtx_300x300_3155nnz

NTIMES=3

#################       MPI Flags        ########################

GMRES_NB_NODES=1
ARNOLDI_NB_NODES=1
MPI_NODES=${shell echo ${GMRES_NB_NODES}+${ARNOLDI_NB_NODES} + 2 | bc}

let "num=$MPI_NODES%24"
let "node=$MPI_NODES/24"

if [ $num -eq 0 ] ; then
  let "NODES=$node"
else
  let "NODES=$node+1"
fi

LSA_GMRES=-lsa_gmres
LSA_ARNOLDI=-lsa_arnoldi

#classical GMRES MPI flags
#MPI_NODES = 2
#################       GMRES Flags       ########################

RESTART_MAX=150
GMRES_PRECISION=1e-10
GMRES_RESTART=${RESTART_MAX}
GMRES_MONITOR=-ksp_monitor_true_residual
KSP_MAX_ITS=20000
PC_TYPE=none
LOG_VIEW=-log_view

#PC_TYPE=jacobi -pc_jacobi_type rowmax -pc_jacobi_abs
#PC_TYPE=asm -pc_asm_block 2 -pc_asm_type none
#PC_TYPE=sor -pc_sor_its 5 -pc_sor_omega 0.179
#CUDA_TYPE=-mat_type aijcusparse -vec_type cuda
#GMRES_FT=-GMRES_FT
GMRES_FLAGS=-ksp_rtol 1e-100 -ksp_divtol 1e1000 -ksp_max_it ${KSP_MAX_ITS} -pc_type ${PC_TYPE} -ksp_atol ${GMRES_PRECISION} -ksp_gmres_restart ${GMRES_RESTART}\
		${GMRES_MONITOR} ${LSA_GMRES} ${GMRES_NB_NODES} -ntimes ${NTIMES} ${CUDA_TYPE}\
		${GMRES_FT}
#		${LOG_VIEW}

#################       ERAM Flags         ########################

ARNOLDI_PRECISION=1e-1
ARNOLDI_NBEIGEN=3
#ARNOLDI_MONITOR=-eps_monitor_conv
#ARNOLDI_FT_SIM=-ArnoldiFT 4

ARNOLDI_FLAGS=-eps_ncv 60 -eps_type arnoldi -eps_true_residual -eps_largest_imaginary -eps_nev ${ARNOLDI_NBEIGEN} -eps_tol ${ARNOLDI_PRECISION} \
                ${ARNOLDI_MONITOR} ${LSA_ARNOLDI} ${ARNOLDI_NB_NODES} -eps_max_it 50 \
		${ARNOLDI_FT_SIM} 
#################       LS Flags           ########################

LS_POWER=10
LS_POLY_APPL=10
LS_LATENCY=1
LS_PC_USE=1
LS_NO_USE_LS=-ksp_ls_nols
LS_HANG_IT=20000
LS_HANG_TIME=100000

#################    LS Reusability Flags  ########################

#LS_CEXPORT=-ksp_ls_cexport
#LS_LOAD=-ksp_ls_load
#LS_LOAD_FILE=./lsqr.bin

LS_FLAGS=-ksp_ls_power ${LS_POWER} ${LS_NO_USE_LS}-ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY} ${LS_CEXPORT} ${LS_LOAD} ${LS_LOAD_FILE}

#################      Flags Combinations      #####################

GLSA_FLAGS=${GMRES_FLAGS} ${ARNOLDI_FLAGS} ${LS_FLAGS} ${DEBUG_KSP_VIEW} ${RUN_FLAGS}

##################################################################
##################     Experimentation     #######################
##################################################################

 yhrun -n ${MPI_NODES} -N ${NODES} ${EXEC} ${GLSA_FLAGS} -mfile ${MDIR}/${MAT}

