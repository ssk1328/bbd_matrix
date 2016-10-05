#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include "bbdf.h"

// gcc t2.c -llapack -std=c99

// Two step process to compute inverse
extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);
// AX=B
extern void dgesv_	(int * n, int * NRHS,double * A, int * LDA, int * IPIV, double * B, int * LDB, int * INFO );
// Matrix Multiplication
//extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
extern void dgemm_ (char * transa, char * transb, int * m, int * n, int * k,
              double * alpha, double * A, int * lda,
              double * B, int * ldb, double * beta,
              double *, int * ldc);
// Scalar Maultiplication
 extern void dscal_ (int * N, double * DA, double * DX, int * INCX);



int print_vector(int size, double * M)
{
	for (size_t i=0; i<size; i++){
		printf(" %g, ", M[i]);
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
	int pivotArray[n];
	int errorHandler;
	double lapackWorkspace[elements]; // change to malloc

	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;

	// matrix inversion M = inv(M)
	dgetrf_(&n, &n, M, &n, pivotArray, &errorHandler);
	dgetri_(&n, M, &n, pivotArray, lapackWorkspace, &elements, &errorHandler);

	// matrix multiplication R = M(nxn)*N(nx1)
	dgemm_(&TRANS, &TRANS, &n, &one, &n, &ALPHA, M, &n, N, &n, &BETA, R, &n);

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
	int pivotArray[m];
	int errorHandler;
	// TODO: Change this to malloc, not a good idea when elements is a large integer
		double lapackWorkspace0[elements];

//	double * lapackWorkspace0;
//	lapackWorkspace0  = (double *) malloc (elements * sizeof(double));

	// matrix inversion
	// Inverted matrix in again stored in A
	dgetrf_(&m, &m, A, &m, pivotArray, &errorHandler);
	dgetri_(&m, A, &m, pivotArray, lapackWorkspace0, &elements, &errorHandler);	

	int one = 1;
	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	double lapackWorkspace1[n*m];
//	double * lapackWorkspace1;

	// matrix multiplication lapackworkspace1 = C(n*m)*inv(A)(m*m)
	dgemm_(&TRANS, &TRANS, &n, &m, &m, &ALPHA, C, &n, A, &m, &BETA, lapackWorkspace1, &n);

	// matrix multiplication tempB = lapackworkspace1(nxm)*B(mxn)
	dgemm_(&TRANS, &TRANS, &n, &n, &m, &ALPHA, lapackWorkspace1, &n, B, &m, &BETA, tempB, &n);

	// TODO: Replace this with matrix vector multiplication
	// matrix multiplication tempG = lapackworkspace1(nxm)*G(mx1)
	dgemm_(&TRANS, &TRANS, &n, &one, &m, &ALPHA, lapackWorkspace1, &n, G, &m, &BETA, tempG, &n);

//	free(lapackWorkspace0);

}


void matrix_xi(int m, int n, double * A, double * B, double * G, double * Xn, double * Xi)
{
	// Xi = A*G - A*B*Xn		
//	int elements = size*size;
	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	int one =1;
	// TODO: Change this to malloc 
	double lapackWorkspaceAB[m*n];
	double lapackWorkspaceX[m*1];

	// matrix multiplication Xi = A*G
    dgemm_(&TRANS, &TRANS, &m, &one, &m, &ALPHA, A, &m, G, &m, &BETA, Xi, &m);

	// matrix multiplication lapackworkspaceAB = A*B
    dgemm_(&TRANS, &TRANS, &m, &n, &m, &ALPHA, A, &m, B, &m, &BETA, lapackWorkspaceAB, &m);

	// matrix multiplication lapackworkspaceX = lapackworkspaceAB*Xn
    dgemm_(&TRANS, &TRANS, &m, &one, &n, &ALPHA, lapackWorkspaceAB, &m, Xn, &n, &BETA, lapackWorkspaceX, &m);

    // Xi = Xi - lapackworkspace2
    matrix_add(m, one, Xi, lapackWorkspaceX, -1.0);
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
	double tempG[n*1];
	for(int i =0; i<(n*1); i++) tempG[i] = 0;

	double AddB[n*n];
	for(int i =0; i<(n*n); i++) AddB[i] = 0;
	double tempB[n*n];
	for(int i =0; i<(n*n); i++) tempB[i] = 0;

	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

//		printf("\n #%d ", i);
//		printf("%s\n", "Matrix Before Inversion");
//		print_matrix(matA[i].nrow, matA[i].ncol, matA[i].matval);

		// tempB = C*inv(A)*B
		// tempG = C*inv(A)*G
		matrix_inv_m(m, n, matA[i].matval, matC[i].matval, matB[i].matval, matG[i].matval, tempB, tempG);

		// AddG = AddG + tempG	:This is vector operation
//		#pragma omp critical
		matrix_add(matGN.nrow, matGN.ncol, AddG, tempG, 1.0);
		// AddB = AddB + tempB	:This is matrix operation
		matrix_add(matAN.nrow, matAN.ncol, AddB, tempB, 1.0);

//		printf("#%d ", i);
//		printf("%s\n", "Matrix After Inversion");
//		print_matrix(matG[i].size, tempG);

	}

//	printf("%s\n", "AddG Matrix WorkSpace After Addition");
//	print_matrix(size, AddG);

//	printf("%s\n", "AddB Matrix WorkSpace After Addition");
//	print_matrix(size, AddB);

//	printf("%s\n", "matA[nMat]");
//	print_matrix(size, matA[nMat].matrix);

//	printf("%s\n", "matG[nMat]");
//	print_matrix(size, matG[nMat].matrix);

	matrix_add(matGN.nrow, matGN.ncol, matGN.matval, AddG, -1.0);	// this should be vector addition
	matrix_add(matAN.nrow, matAN.ncol, matAN.matval, AddB, -1.0);

//	printf("%s\n", "After sub matA[nMat]");
//	print_matrix(size, matA[nMat].matrix);

//	printf("%s\n", "After sub matG[nMat]");
//	print_matrix(size, matG[nMat].matrix);

	matrix_invm(matAN.nrow, matAN.matval, matGN.matval, matXN.matval);

//	printf("%s\n", "After final Inversion matX[nMat]");
//	print_matrix(size, matX[nMat].matrix);

	for (int i=0; i<nMat; i++){

		matrix_xi(m, n, matA[i].matval, matB[i].matval, matG[i].matval, matXN.matval, matX[i].matval);
//		printf("%s\n", "-------------------------------------------------" );
//		printf("%d th ",i);
//		printf("%s\n", "Solution Xi Matrices");
		print_matrix(matX[i].nrow, matX[i].ncol,  matX[i].matval);

	}

}
