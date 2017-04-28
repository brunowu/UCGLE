/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "Libs/mpi_lsa_com.h"
#include "slepceps.h"
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include "./Solvers/GMRES/gmres.h"
#include "./Solvers/Arnoldi/arnoldi.h"
#include "./Libs/read_matrix.h"
#include "./Solvers/LS/lsqr.h"
#include "./Solvers/Father/father.h"
#include "./Libs/real2complex.h"
