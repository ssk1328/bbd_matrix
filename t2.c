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

int main() {
	int N = 4;
	double M[16] = { 1, 0, 0, 0,
					 0, 1, 0, 0,
                   	 0, 0, 1, 0,
                   	 0, 0, 0, 1 };

//    int pivotArray[4]; //since our matrix has three rows
//    int errorHandler;
//    double lapackWorkspace[16];

//    dgetrf_(&N, &N, M, &N, pivotArray, &errorHandler);
//    dgetri_(&N, M, &N, pivotArray, lapackWorkspace, &NN, &errorHandler);

	matrix_inv(N, M);

    for (size_t row = 0; row < N; ++row)
    {   for (size_t col = 0; col < N; ++col)
        {   printf ("%g", M[row*N+col]);
            if (N-1 != col)
            {   printf (", ");   }   }
        if (N-1 != row)
        {   printf ("\n");   }   }
    return 0;   
}

