#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include "bbdf.h"
#include "lapacke.h"
#include "cblas.h"

// gcc t2.c -llapack -std=c99

// Two step process to compute inverse
//extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
//extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);
// AX=B
//extern void dgesv_	(int * n, int * NRHS,double * A, int * LDA, int * IPIV, double * B, int * LDB, int * INFO );
// Matrix Multiplication
//extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
extern void dgemm_ (char * transa, char * transb, int * m, int * n, int * k, double * alpha, double * A, int * lda, double * B, int * ldb, double * beta, double *, int * ldc);
// Scalar Maultiplication
// extern void dscal_ (int * N, double * DA, double * DX, int * INCX);

int print_matrix(int nrow, int ncol, double * M) {

	for (int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			printf("%.3f ", M[i*ncol + j]);
		}
		printf("\n");
	}
	return 0;
}

void matrix_add(int row, int col, double * A, double * B, double scalar)
{
	// A = (A)+(s*B)
	int elements = row*col;
	for (size_t i=0; i<elements; i++){
		A[i] = A[i] + scalar * B[i];
	}
}

void vector_add(int row, int col, double * A, double * B, double scalar)
{	
	// A = (A)+(s*B)
	int N = row*col;
	for (size_t i=0; i<N; i++){
		A[i] = A[i] + scalar * B[i];
	}
}

void matrix_invm (int n, double * M, double * N, double * R) 
{
	//	Agruments:
	//	n: Size of matrix M is nxn 
	// 	M: Pointer to start of Matrix, stored in an array

	// R = inv(M)*N
	int one = 1 ;
	int elements = n*n;
	int * ipiv = (int*) malloc(n*sizeof(int));
	int info;

	char TRANSN = 'N';
	char TRANST = 'T';
	double ALPHA = 1.0;
	double BETA = 0.0;

	// matrix inversion M = inv(M)
	//	dgetrf_(&n, &n, M, &n, pivotArray, &errorHandler);
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, M, n, ipiv);
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, M, n, ipiv); //, lapackWorkspace0, elements);

//	C(rowsA x colsB) = A(rowsA x common)*B(common x colsB)
//	dgemm_(&TRANSN, &TRANSN, &colsB, &rowsA, &common, 	&ALPHA, B, &colsB, A, &common, 	&BETA, C, &colsB);
//	matrix multiplication R = M(nxn)*N(nx1)
	dgemm_(&TRANSN, &TRANSN, &one, &n, &n, 				&ALPHA, N, &one, M, &n, 		&BETA, R, &one)	;
}


void matrix_inv_m (int m, int n, double * A, double * C, double * B, double * G, double * tempB, double * tempG) {

	//	Agruments:
	//  m:		int - A is mxm matrix
	//  n:		int - B is mxn matrix
	//			int - C is nxm matrix
	// 	A: 		pointer to double - Pointer to start of Matrix A, stored in contagious array
	// 	C: 		pointer to double - Pointer to start of Matrix B, stored in contagious array
	// 	B: 		pointer to double - Pointer to start of Matrix C, stored in contagious array
	// 	G: 		pointer to double - Pointer to start of Matrix D, stored in contagious array

	// tempB  = C*inv(A)*B
	// tempG  = C*inv(A)*G

	int elements = m*m;
	int * ipiv = (int*) malloc(m*sizeof(int));
	int info;
	double * lapackWorkspace0 = (double *) malloc (elements * sizeof(double));

	// matrix inversion
	// Inverted matrix in again stored in A
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, A, m, ipiv);
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, m, A, m, ipiv);

	int one = 1;
	char TRANSN = 'N';
	char TRANST = 'T';
	double ALPHA = 1.0;
	double BETA = 0.0;
	double * lapackWorkspace1 = (double *) malloc (n*m * sizeof(double));

	// matrix multiplication lapackworkspace1 = C(n*m)*inv(A)(m*m)
//	dgemm_(&TRANSN, &TRANST, &n, &m, &m, &ALPHA, C, &n, A, &m, &BETA, lapackWorkspace1, &n);
//	C(rowsA x colsB) = A(rowsA x common)*B(common x colsB)
//	dgemm_(&TRANSN, &TRANSN, &colsB, &rowsA, &common, &ALPHA, B, &colsB, A, &common, &BETA, C, &colsB);
	dgemm_(&TRANSN, &TRANSN, &m, &n, &m, &ALPHA, A, &m, C, &m, &BETA, lapackWorkspace1, &m);

	// matrix multiplication tempB = lapackworkspace1(nxm)*B(mxn)
	dgemm_(&TRANSN, &TRANSN, &n, &n, &m, &ALPHA, B, &n, lapackWorkspace1, &m, &BETA, tempB, &n);

	// matrix multiplication tempG = lapackworkspace1(nxm)*G(mx1)
	dgemm_(&TRANSN, &TRANSN, &one, &n, &m, &ALPHA, G, &one, lapackWorkspace1, &m, &BETA, tempG, &one);

	free(lapackWorkspace0);
	free(lapackWorkspace1);
}

void matrix_xi(int m, int n, double * A, double * B, double * G, double * Xn, double * Xi)
{
//	Xi = A*G - A*B*Xn
	char TRANSN = 'N';
	char TRANST = 'T';
	double ALPHA = 1.0;
	double BETA = 0.0;
	int one = 1;
	double * lapackWorkspaceAB = (double*) malloc ((m*n)*sizeof(double));
	double * lapackWorkspaceX = (double*) malloc ((m*1)*sizeof(double));

//	matrix multiplication Xi = A(mxm)*G(mx1)
//	dgemm_(&TRANST, &TRANSN, &m, &one, &m, &ALPHA, A, &m, G, &m, &BETA, Xi, &m);
//	C(rowsA x colsB) = A(rowsA x common)*B(common x colsB)
//	dgemm_(&TRANSN, &TRANSN, &colsB, &rowsA, &common, &ALPHA, B, &colsB, A, &common, &BETA, C, &colsB);
	dgemm_(&TRANSN, &TRANSN, &one, &m, &m, &ALPHA, G, &one, A, &m, &BETA, Xi, &one);

//	matrix multiplication lapackworkspaceAB(mxn) = A(mxm)*B(mxn)
//	dgemm_(&TRANSN, &TRANSN, &m, &n, &m, &ALPHA, A, &m, B, &m, &BETA, lapackWorkspaceAB, &m);
	dgemm_(&TRANSN, &TRANSN, &n, &m, &m, &ALPHA, B, &n, A, &m, &BETA, lapackWorkspaceAB, &n);

// 	matrix multiplication lapackworkspaceX = lapackworkspaceAB(mxn)*Xn(nx1)
//	dgemm_(&TRANSN, &TRANSN, &m, &one, &n, &ALPHA, lapackWorkspaceAB, &m, Xn, &n, &BETA, lapackWorkspaceX, &m);
	dgemm_(&TRANSN, &TRANSN, &one, &m, &n, &ALPHA, Xn, &one, lapackWorkspaceAB, &n, &BETA, lapackWorkspaceX, &one);

	// Xi = Xi - lapackworkspace2
	matrix_add(m, one, Xi, lapackWorkspaceX, -1.0);

	free(lapackWorkspaceAB);
	free(lapackWorkspaceX);
}


void solve_bbd( int nMat, int m, int n,
				struct Matrix * matA,  struct Matrix * matB, 
				struct Matrix * matC,  struct Matrix * matG, 
				struct Matrix * matX,  struct Matrix matAN, 
				struct Matrix   matGN, struct Matrix matXN)
{

/*
	A 		   B    X     G
	   A 	   B    X     G
	      A    B  * X  =  G
		     A B    X     G
	C  C  C  C AN   XN    GN

	Arguments:
	nMat: 	Number of non border diagonal blocks (nMat = 4 in above example)
	matA:	Pointer to array of Struct Matrix A, diagonal matrix blocks
	matB:	Pointer to array of Struct Matrix B, bottom border matrix blocks
	matC:	Pointer to array of Struct Matrix C, right border matrix blocks
	matG:	Pointer to array of Struct Matrix G, vector
	matX:	Pointer to array of Struct Matrix X, vector
*/

	double AddG[n*1];
	for(int i =0; i<(n*1); i++) AddG[i] = 0;
	double tempG[nMat][n*1];

	double AddB[n*n];
	for(int i =0; i<(n*n); i++) AddB[i] = 0;
	double tempB[nMat][n*n];

    clock_t tstart = clock();

//	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		// tempB = C*inv(A)*B
		// tempG = C*inv(A)*G
		matrix_inv_m(m, n, matA[i].matval, matC[i].matval, matB[i].matval, matG[i].matval, tempB[i], tempG[i]);
	}

//	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		// AddG = AddG + tempG	:This is vector operation
		matrix_add(matGN.nrow, matGN.ncol, AddG, tempG[i], 1.0);
		// AddB = AddB + tempB	:This is matrix operation
		matrix_add(matAN.nrow, matAN.ncol, AddB, tempB[i], 1.0);
	}

    clock_t p1 = clock();

    printf("Part1 Time taken: %.4fms\n", (double)(p1 - tstart)/(CLOCKS_PER_SEC*0.001));

	matrix_add(matGN.nrow, matGN.ncol, matGN.matval, AddG, -1.0);
	matrix_add(matAN.nrow, matAN.ncol, matAN.matval, AddB, -1.0);

	matrix_invm(matAN.nrow, matAN.matval, matGN.matval, matXN.matval);

    clock_t p2 = clock();
    printf("Part2 Time taken: %.4fms\n", (double)(p2 - p1)/(CLOCKS_PER_SEC*0.001));

	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		matrix_xi(m, n, matA[i].matval, matB[i].matval, matG[i].matval, matXN.matval, matX[i].matval);
//		printf("%s\n", "-------------------------------------------------" );
//		printf("%d th ",i);
//		printf("%s\n", "Solution Xi Matrices");
//		print_matrix(matX[i].nrow, matX[i].ncol,  matX[i].matval);
	}
//	print_matrix(matXN.nrow, matXN.ncol,  matXN.matval);

    clock_t p3 = clock();
    printf("Part3 Time taken: %.4fms\n", (double)(p3 - p2)/(CLOCKS_PER_SEC*0.001));
}

void solve_bbd_full (int N, int m, int n,
				struct Matrix * matA, struct Matrix * matB, 
				struct Matrix * matC, struct Matrix * matG, 
				struct Matrix  matAN, struct Matrix  matGN){

	int size = (N-1)*m+n;
	double * A = (double*) malloc(size*size*sizeof(double));
	double * AT = (double*) malloc(size*size*sizeof(double));
	double * G = (double*) malloc(size*sizeof(double));
	double * X = (double*) malloc(size*sizeof(double));
	int i,j,k;

	int nrhs = 1;
	int lda = size;
	int ldb = nrhs;
	int * ipiv = (int*) malloc(size*sizeof(int));
	int info;

	for(i=0; i<size*size; i++){
		A[i]=0;
	}

//	printf("Size is :%d\n", size );

	for(i = 0; i<N-1; i++){
		for(j=0; j<m; j++){
			for(k=0; k<m; k++){
				A[ (i*m+j)*size + (i*m)+k] = matA[i].matval[j*m+k];
			}
		}

		for(j=0; j<n; j++){
			for(k=0; k<m; k++){
				A[ ((N-1)*m+j)*size+ (i*m)+k] = matC[i].matval[j*m+k];
			}
		}

		for(j=0; j<m; j++){
			for(k=0; k<n; k++){
				A[(i*m+j)*size + (N-1)*m+k] = matB[i].matval[j*n+k];
			}
		}

		for(j=0; j<m; j++){
			G[i*m+j] = matG[i].matval[j];
		}
	}

	for(j=0; j<n; j++){
		for(k=0; k<n; k++){
			A[((N-1)*m+j)*size + (N-1)*m+k] = matAN.matval[j*n+k];
		}
	}

	for(j=0; j<n; j++){
			G[(N-1)*m+j]= matGN.matval[j];
	}

//	printf("--------------------------------------------------------------------\n" );
//	printf("Full matrix A: \n" );
//	print_matrix(size, size, A);
//	printf("Full matrix G: \n" );
//	print_matrix(size, 1, G);

	clock_t t1 = clock();
//	dgesv_(&size, &one, AT, &size, pv, G, &size, &INFO);
	info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, size, nrhs, A, lda, ipiv, G, ldb);
	clock_t t2 = clock();

	printf("INFO VALUE: %d\n", info);

	printf("Total Time taken for full Matrix: %.4fms\n", (double)(t2 - t1)/(CLOCKS_PER_SEC*0.001));
//	printf("Print Xi solutions from full matrix \n");
//	print_matrix(size, 1, G);

}