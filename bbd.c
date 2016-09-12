#include <stdio.h>
#include <stddef.h>

// gcc t2.c -llapack -std=c99

extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);

extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);

void matrix_inv (int size, double * M, double * B) {

	//	Agruments:
	//	size: Size of matrix
	// 	M: Pointer to start of Matrix, stored in an array

	// Inverted matrix in again stored in M

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	double lapackWorkspace[elements];

	dgetrf_(&size, &size, M, &size, pivotArray, &errorHandler);
	dgetri_(&size, M, &size, pivotArray, lapackWorkspace, &elements, &errorHandler);	
	
	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;

    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, A, &size, B, &size, &BETA, B, &size);

}

int print_matrix(int size, double * M)
{
	/* code */
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

struct Matrix {
	int size;
	double matrix[100];
	// row first form 
	// Each matrix can have max size of 10x10
};


void solve_bbd( int nMat, struct Matrix * mat){

	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		printf("\n #%d ", i);
		printf("%s\n", "Matrix Before Inversion");
		print_matrix(mat[i].size, mat[i].matrix);

		matrix_inv(mat[i].size, matA[i].matrix, matB[i].matrix);

		printf("#%d ", i);	
		printf("%s\n", "Matrix After Inversion");
		print_matrix(mat[i].size, mat[i].matrix);
	}

}

int main() {

	int nMat = 6;
	struct Matrix matA[nMat];
	struct Matrix matB[nMat-1];
	struct Matrix matC[nMat-1];
	struct Matrix matG[nMat];

	int n; // temp

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
	}

	solve_bbd(nMat, mat);

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
	printf("\n %s\n", "End End End End End End End End End End End End" );
    return 0;   
}

