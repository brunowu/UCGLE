/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */

#include "mat_generate.h"
#include <interface/C/c_wrapper.h>

PetscErrorCode mat_generate(Mat * A, MPI_Comm comm){

	int *rows, *cols;
        int size_row, size_col;
        double *real, *imag;

	char spectra[PETSC_MAX_PATH_LEN];
	PetscBool flgsmg2s, flglb, flgsize, flgno, flgsptr;
	PetscErrorCode ierr;
	
	PetscInt lbandwidth, size, nbOne;
	ierr= PetscOptionsHasName(NULL,NULL,"-smg2s",&flgsmg2s);CHKERRQ(ierr);
	ierr= PetscOptionsGetInt(NULL,NULL,"-lbandwidth",&lbandwidth,&flglb);CHKERRQ(ierr);
        ierr= PetscOptionsGetInt(NULL,NULL,"-size",&size,&flgsize);CHKERRQ(ierr);
        ierr= PetscOptionsGetInt(NULL,NULL,"-continous1",&nbOne,&flgno);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-sptr",spectra,PETSC_MAX_PATH_LEN-1,&flgsptr);CHKERRQ(ierr);
	if(flgsmg2s) PetscPrintf(comm, "Info ]> Using SMG2S to generate test matrices\n");
	if(!flglb){
		lbandwidth = 3;
		PetscPrintf(comm, "Remainder ]> low bandwidth in SMG2S is not set, use the default value = 3\n");
	}

        if(!flgsize){
                size = 100;
                PetscPrintf(comm, "Remainder ]> matrix size in SMG2S is not set, use the default value = 100\n");
        }

        if(!flgno){
                nbOne = 2;
                PetscPrintf(comm, "Remainder ]> Continus 1 in SMG2S is not set, use the default value = 2\n");
        }

	if(!flgsptr){
                strncpy(spectra, " ", sizeof(spectra));
		PetscPrintf(comm, "Remainder ]> given spectra file in SMG2S is not set, use the default method to generate\n");
        }

	struct NilpotencyInt *n;
        n = newNilpotencyInt();
        NilpType1(n, nbOne, size);
        showNilpotencyInt(n);

        struct parMatrixSparseComplexDoubleInt *m;
        m = newParMatrixSparseComplexDoubleInt();

	smg2sComplexDoubleInt(m, size, n, lbandwidth ,spectra, comm);

        Loc_ConvertToCSRComplexDoubleInt(m);
        Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);

        rows = (int *) malloc(size_row*sizeof(int));
        cols = (int *) malloc(size_col*sizeof(int));

        real = (double *) malloc(size_col*sizeof(double));
        imag = (double *) malloc(size_col*sizeof(double));

        Loc_CSRGetRowsArrays(m, size_row, &rows, size_col, &cols, &real, &imag);

        PetscInt *i, *j, q, p;
        PetscScalar *a;

        PetscMalloc1((PetscInt)size_row, &i);
        PetscMalloc1((PetscInt)size_col, &j);
        PetscMalloc1((PetscInt)size_col, &a);


	for(q =0; q < size_row; q++){
                i[q] = rows[q];
        }

        for(p = 0; p < size_col; p++){
                j[p] = cols[p];
                a[p] = (PetscReal) real[p] + PETSC_i * (PetscReal) imag[p];
	}


        MatCreate(PETSC_COMM_WORLD, A);
        MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, size_row-1, size_row-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, A);


//	free(rows); free(cols); free(real); free(imag);
//	PetscFree(i); PetscFree(j); PetscFree(a);

//	ReleaseParMatrixSparseComplexDoubleInt(&m);
	ReleaseNilpotencyInt(&n);
	
	return ierr;

}

