#include <stdio.h>
#include <stddef.h>

// gcc t2.c -llapack -std=c99

// lapack and blas functions used

// Two step process to compute inverse
extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);

int print_matrix(int size, double * M)
{
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

void matrix_inv_m (int size, double * A) {

	//	Agruments:
	//	size: 	int - Size of matrix
	// 	A: 		pointer to double - Pointer to start of Matrix M, stored in contagious array

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	// TODO: Change this to malloc, not a good idea when elements is a large integer
	double lapackWorkspace0[elements];

	// matrix inversion
	// Inverted matrix in again stored in A
	dgetrf_(&size, &size, A, &size, pivotArray, &errorHandler);
	dgetri_(&size, A, &size, pivotArray, lapackWorkspace0, &elements, &errorHandler);	
}


struct Matrix {
	// row first form 
	// Each matrix can have max size of 10x10
	// TODO Add support for having non square matrices
	int nrow;
	int ncol;
	double matval[10][10];
};

int main() {

	int N = 3; 	// number of blocks in diagonal
	int m = 4;	// mxm is the size of Ai
	int n = 3;	// nxn is the size of AN
	
	struct Matrix matA[N-1];
	struct Matrix matB[N-1];
	struct Matrix matC[N-1];
	struct Matrix matG[N-1];
	struct Matrix matAN;
	struct Matrix matGN;

	for (int i=0; i<N-1; i++){

		matA[i].nrow = m;
		matA[i].ncol = m;
		matB[i].nrow = m;
		matB[i].ncol = n;
		matC[i].nrow = n;
		matC[i].ncol = m;
	}

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


	return 0;
}

