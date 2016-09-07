#include <stdio.h>
#include <stddef.h>

// gcc t2.c -llapack -std=c99

extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);


void matrix_inv (int size, double * M) {

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

}

int print_matrix(int size, double * M)
{
	/* code */
	int N = size;
	int elements = size*size;

	for (size_t i=0; i<elements; i++){
		printf("%g, ", M[i]);
		if( i%N == (N-1)) {
			printf("\n");
		}
	}

	return 0;
}


int main() {
	int N = 4;
	double M[16] = { 1, 0, 0, 0,
					 0, 1, 0, 0,
                   	 0, 0, 1, 0,
                   	 0, 0, 0, 1 };

	matrix_inv(N, M);
	print_matrix(N,M);

    return 0;   
}

