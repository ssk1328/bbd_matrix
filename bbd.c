#include <stdio.h>
#include <stddef.h>

// gcc t2.c -llapack -std=c99

// lapack and blas functions used

// Two step process to compute inverse
extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);
// AX=B
extern void dgesv_	(int * n, int * NRHS,double * A, int * LDA, int * IPIV, double * B, int * LDB, int * INFO )
// Matrix Multiplication
extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
// Scalar Maultiplication
extern void dscal_ (int * N, double * DA, double * DX, int * INCX);

int print_matrix(int size, double * M) {
	int N = size;
	int elements = size*size;
	for (size_t i=0; i<elements; i++){
		printf(" %g, ", M[i]);
		if( i%N == (N-1)) {
			printf("\n");
		}
	}
	return 0;
}

int print_vector(int size, double * M)
{
	for (size_t i=0; i<size; i++){
		printf(" %g, ", M[i]);
		printf("\n");		
	}
	return 0;
}

void matrix_add(int size, double * A, double * B, double scalar)
{
	// A = (A)+(s*B)
	int N = size;
	int elements = size*size;
	for (size_t i=0; i<elements; i++){
		A[i] = A[i] + scalar * B[i];
	}
}

void vector_add(int size, double * A, double * B, double scalar)
{	
	// A = (A)+(s*B)
	int N = size;
	for (size_t i=0; i<N; i++){
		A[i] = A[i] + scalar * B[i];
	}
}

void matrix_invm (int size, double * M, double * N, double * R) 
{
	//	Agruments:
	//	size: Size of matrix
	// 	M: Pointer to start of Matrix, stored in an array

	// R = inv(M)*N

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	double lapackWorkspace[elements]; // change to malloc

	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;

	// matrix inversion M = inv(M)
	dgetrf_(&size, &size, M, &size, pivotArray, &errorHandler);
	dgetri_(&size, M, &size, pivotArray, lapackWorkspace, &elements, &errorHandler);	

	// matrix multiplication R = M*N
	dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, M, &size, N, &size, &BETA, R, &size);

}

void matrix_inv_m (int size, double * A, double * C, double * B, double * G, double * tempB, double * tempG) {

	//	Agruments:
	//	size: 	int - Size of matrix
	// 	A: 		pointer to double - Pointer to start of Matrix M, stored in contagious array
	// 	C: 		pointer to double - Pointer to start of Matrix B, stored in contagious array
	// 	B: 		pointer to double - Pointer to start of Matrix C, stored in contagious array
	// 	G: 		pointer to double - Pointer to start of Matrix D, stored in contagious array

	// tempB  = C*inv(A)*B
	// tempG  = C*inv(A)*G

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	// TODO: Change this to malloc, not a good idea when elements is a large integer
	double lapackWorkspace0[elements];

	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;

	// matrix inversion
	// Inverted matrix in again stored in A
	dgetrf_(&size, &size, A, &size, pivotArray, &errorHandler);
	dgetri_(&size, A, &size, pivotArray, lapackWorkspace0, &elements, &errorHandler);	

	// matrix multiplication lapackworkspace0 = C*inv(A)
	dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, C, &size, A, &size, &BETA, lapackWorkspace0, &size);

	// matrix multiplication tempB = lapackworkspace0*B
	dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, lapackWorkspace0, &size, B, &size, &BETA, tempB, &size);

	// TODO: Replace this with matrix vector multiplication
	// matrix multiplication tempG = lapackworkspace0*G
	dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, lapackWorkspace0, &size, G, &size, &BETA, tempG, &size);

}


void matrix_xi(int size, double * A, double * B, double * G, double * Xn, double * Xi)
{
	// Xi = A*G - A*B*Xn		
	int elements = size*size;
	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	// TODO: Change this to malloc 
	double lapackWorkspace1[elements];
	double lapackWorkspace2[elements];

	// matrix multiplication Xi = A*G
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, A, &size, G, &size, &BETA, Xi, &size);

	// matrix multiplication lapackworkspace1 = A*B
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, A, &size, B, &size, &BETA, lapackWorkspace1, &size);

	// matrix multiplication lapackworkspace2 = lapackworkspace1*Xn
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, lapackWorkspace1, &size, Xn, &size, &BETA, lapackWorkspace2, &size);

    // Xi = Xi - lapackworkspace2
    matrix_add(size, Xi, lapackWorkspace2, -1.0);
}

struct Matrix {
	// row first form 
	// Each matrix can have max size of 10x10
	// TODO Add support for having non square matrices
	int size;
	double matrix[100];
};

void solve_bbd( int nMat, struct Matrix * matA, struct Matrix * matB, struct Matrix * matC, struct 	Matrix * matG, struct 	Matrix * matX) {
/*

	A 		   B   X   G
	   A 	   B   X   G
	      A    B * X = G
		     A B   X   G
	C  C  C  C A   X   G

	Arguments:
	nMat: 	Number of diagonal blocks minus one (nMat = 4 in above example)
	matA:	Pointer to array of Struct Matrix A, diagonal matrix blocks
	matB:	Pointer to array of Struct Matrix B, bottom border matrix blocks
	matC:	Pointer to array of Struct Matrix C, right border matrix blocks
	matG:	Pointer to array of Struct Matrix G, vector
	matX:	Pointer to array of Struct Matrix X, vector

*/

	int size = matA[0].size;
	int elements = size*size;

	double AddG[elements];
	for(int i =0; i<elements; i++) AddG[i] = 0;
	double tempG[elements];
	for(int i =0; i<elements; i++) tempG[i] = 0;

	double AddB[elements];
	for(int i =0; i<elements; i++) AddB[i] = 0;
	double tempB[elements];
	for(int i =0; i<elements; i++) tempB[i] = 0;

	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		printf("\n #%d ", i);
		printf("%s\n", "Matrix Before Inversion");
		print_matrix(matA[i].size, matA[i].matrix);

		// tempB = C*inv(A)*B
		// tempG = C*inv(A)*G
		matrix_inv_m(matA[i].size, matA[i].matrix, matC[i].matrix, matB[i].matrix, matG[i].matrix, tempB, tempG);

		// AddG = AddG + tempG	:This is vector operation
		matrix_add(matA[i].size, AddG, tempG, 1.0);
		// AddB = AddB + tempB	:This is matrix operation
		matrix_add(matA[i].size, AddB, tempB, 1.0);

		printf("#%d ", i);
		printf("%s\n", "Matrix After Inversion");
		print_matrix(matG[i].size, tempG);

	}

//	printf("%s\n", "AddG Matrix WorkSpace After Addition");
//	print_matrix(size, AddG);

	printf("%s\n", "AddB Matrix WorkSpace After Addition");
	print_matrix(size, AddB);

	printf("%s\n", "matA[nMat]");
	print_matrix(size, matA[nMat].matrix);

	printf("%s\n", "matG[nMat]");
	print_matrix(size, matG[nMat].matrix);

	matrix_add(size, matG[nMat].matrix, AddG, -1.0);	// this should be vector addition
	matrix_add(size, matA[nMat].matrix, AddB, -1.0);

	printf("%s\n", "After sub matA[nMat]");
	print_matrix(size, matA[nMat].matrix);

	printf("%s\n", "After sub matG[nMat]");
	print_matrix(size, matG[nMat].matrix);

	matrix_invm(size, matA[nMat].matrix, matG[nMat].matrix, matX[nMat].matrix);

	printf("%s\n", "After final Inversion matX[nMat]");
	print_matrix(size, matX[nMat].matrix);

	for (int i=0; i<nMat; i++){

		matrix_xi(size, matA[i].matrix, matB[i].matrix, matG[i].matrix, matX[nMat].matrix, matX[i].matrix);
		printf("%s\n", "-------------------------------------------------" );
		printf("%d th ",i);
		printf("%s\n", "Solution Xi Matrices");
		print_matrix(size, matX[i].matrix);

	}

}

int main() {

	int nMat = 2;
	struct Matrix matA[nMat+1];
	struct Matrix matB[nMat];
	struct Matrix matC[nMat];
	struct Matrix matG[nMat+1];
	struct Matrix matX[nMat+1];

	int n; // temp

	for (int i=0; i<nMat+1; i++){

		matA[i].size = 3;
		matB[i].size = 3;
		matC[i].size = 3;
		matG[i].size = 3;
		n = matA[i].size;

		// Initialize each Matrix

		if (i<nMat)
		{
			/* code */
			for (int j=0; j<n*n; j++) {
				if( j/n == j%n ) {
					matA[i].matrix[j] = (i+1)*2;
					matB[i].matrix[j] = (i+1)*2;
					matC[i].matrix[j] = (i+1)*2;
					matG[i].matrix[j] = (i+1)*2;
				}
				else {
					matA[i].matrix[j] = 0;
					matB[i].matrix[j] = 0;
					matC[i].matrix[j] = 0;
					matG[i].matrix[j] = 0;
				}
			}
		}
		else
		{
			for (int j=0; j<n*n; j++) {
					if( j/n == j%n ) {
						matA[i].matrix[j] = (i+1)*4;
						matG[i].matrix[j] = (i+1)*4;
					}
					else {
						matA[i].matrix[j] = 0;
						matG[i].matrix[j] = 0;
					}
				}	
		}
	}

	solve_bbd(nMat, matA, matB, matC, matG, matX);

/*
	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		mat[i].size = i+3;
		n = mat[i].size;

		// Initialize each Matrix
		for (int j=0; j<n*n; j++) {
			if( j/n == j%n )
				mat[i].matrix[j] = (i+1)*2;
			else 
				mat[i].matrix[j] = 0;
		}

	printf("\n #%d ", i);
	printf("%s\n", "Matrix Before Inversion");
	print_matrix(mat[i].size, mat[i].matrix);
	matrix_inv(mat[i].size, mat[i].matrix);
	printf("#%d ", i);	
	printf("%s\n", "Matrix After Inversion");
	print_matrix(mat[i].size, mat[i].matrix);

	}
*/
	
/*
	int N1 = 4;
	double M1[16] = { 1, 0, 0, 0,
					  0, 1, 0, 0,
                   	  0, 0, 1, 0,
                   	  0, 0, 0, 1 
                   	};

	int N2 = 4;
	double M2[16] = {  1, 0, 0, 
					  0, 1, 0, 
                   	  0, 0, 1 
                   };


	matrix_inv(N2, M2);
	print_matrix(N2,M2);
*/
	printf("\n %s\n", "-------------------------------------------------" );
	return 0;
}

