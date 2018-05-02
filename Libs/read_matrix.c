/*Copyright (c) 2011—2017. Pierre-Yves AQUILANTI and Xinzhe WU in Maison de la Simulation. All rights reserved */
#include "read_matrix.h"
#include <interface/C/c_wrapper.h>
PetscErrorCode read_matrix_vector(Mat * A, Vec * v, int *comm){
	char filea[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flaga,flagb, flgsmg2s;
	PetscViewer fd;
	PetscInt size,sizea;
	PetscScalar scal;
	PetscReal vnorm;

	ierr= PetscOptionsHasName(NULL,NULL,"-smg2s",&flgsmg2s);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);
	if (!flaga) {
		if(!flgsmg2s){
		        sprintf(err,"Error : mfile is not properly set -> %s\n",filea);
	                SETERRQ(*comm,(PetscErrorCode) 83,err);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD,"Remainder ]> USING SMG2S to generate test matrix\n");
			mat_generate(A, *comm);
		}
	}
	else{
		/* read matrix file */
	        PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
       		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
        	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
        	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);

	}

	/* read matrix file */
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);
/*
	int rank;
	PetscInt start, end;

	int pp;

	MPI_Comm_rank(*comm, &rank);
	MPI_Comm_size(*comm, &pp);
	MatGetOwnershipRange(*A, &start, &end);
	printf("RANK %d / Size %d: Loaded Matrix of size : %d %d, ownership: start = %d, end = %d\n",rank,pp, size,sizea, start, end);
//	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",size,sizea);

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);

	int rs, cs;
        int *rows, *cols;
        int size_row, size_col;
        double *real, *imag;

	struct NilpotencyInt *n;
        n = newNilpotencyInt();
        NilpType1(n, 2, 10);
        showNilpotencyInt(n);
     
        struct parMatrixSparseComplexDoubleInt *m;
	m = newParMatrixSparseComplexDoubleInt();

        smg2sComplexDoubleInt(m, 10, n, 3 ," ", *comm);

	LOC_MatViewComplexDoubleInt(m);
        Loc_ConvertToCSRComplexDoubleInt(m);
	Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);

        rows = (int *) malloc(size_row*sizeof(int));
        cols = (int *) malloc(size_col*sizeof(int));

        real = (double *) malloc(size_col*sizeof(double));
        imag = (double *) malloc(size_col*sizeof(double));

        Loc_CSRGetRowsArrays(m, size_row, &rows, size_col, &cols, &real, &imag);
        printf("helloooo\n");

        PetscInt *i, *j, q, p;
        PetscScalar *a;

        Mat S;

        PetscMalloc1((PetscInt)size_row, &i);
        PetscMalloc1((PetscInt)size_col, &j);
        PetscMalloc1((PetscInt)size_col, &a);

        for(q =0; q < size_row; q++){
                i[q] = rows[q];
	}

	for(p = 0; p < size_col; p++){
                j[p] = cols[p];
                a[p] = real[p] + PETSC_i * imag[p];
	}


	MatCreate(PETSC_COMM_WORLD, &S);
        MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, size_row-1, size_row-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &S);
	MatView(S, PETSC_VIEWER_STDOUT_WORLD);

	free(rows);free(cols); free(real);free(imag);
	PetscFree(i);PetscFree(j);PetscFree(a);
*/

	if (!flagb) {
		/* the user did not provide a vector, so generate it*/
		generate_random_seed_vector(size, -10.0, 10.0, 0, v);
		VecNorm(*v, NORM_2, &vnorm);
		VecScale(*v, 0.01/vnorm);		
		PetscPrintf(PETSC_COMM_WORLD,"Generated right hand side matrix b\n");		
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",fileb);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileb,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		VecCreate(PETSC_COMM_WORLD,v);
		ierr=VecLoad(*v,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

	return 0;
}
