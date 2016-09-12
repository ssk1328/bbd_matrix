#include <stdio.h>
#include <stddef.h>

// gcc t2.c -llapack -std=c99

extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);
extern void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
extern void dscal_ (int * N, double * DA, double * DX, int * INCX);

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

void matrix_invm (int size, double * M, double * N, double * R) {

	//	Agruments:
	//	size: Size of matrix
	// 	M: Pointer to start of Matrix, stored in an array

	// M = inv(M)*N

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	double lapackWorkspace[elements]; // change to malloc

	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;

	// matrix inversion
	dgetrf_(&size, &size, M, &size, pivotArray, &errorHandler);
	dgetri_(&size, M, &size, pivotArray, lapackWorkspace, &elements, &errorHandler);	

	// matrix multiplication lapackworkspace = B*C
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, M, &size, N, &size, &BETA, R, &size);

}

void matrix_inv_m (int size, double * M, double * B, double * C, double * D) {

	//	Agruments:
	//	size: Size of matrix
	// 	M: Pointer to start of Matrix, stored in an array

	// Inverted matrix in again stored in M

	int elements = size*size;
	int pivotArray[size];
	int errorHandler;
	double lapackWorkspace[elements];
	double lapackWorkspace2[elements];

	char TRANS = 'N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	// matrix inversion
	dgetrf_(&size, &size, M, &size, pivotArray, &errorHandler);
	dgetri_(&size, M, &size, pivotArray, lapackWorkspace, &elements, &errorHandler);	

	// matrix multiplication lapackworkspace = B*C
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, B, &size, C, &size, &BETA, lapackWorkspace, &size);

	// matrix multiplication lapackworkspace2 = B*D
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, B, &size, D, &size, &BETA, lapackWorkspace2, &size);

	// matrix multiplication C = M*lapackworkspace
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, M, &size, lapackWorkspace, &size, &BETA, C, &size);

	// matrix multiplication D = M*lapackworkspace2
    dgemm_(&TRANS, &TRANS, &size, &size, &size, &ALPHA, M, &size, lapackWorkspace2, &size, &BETA, D, &size);

}

void matrix_add(int size, double * A, double * B, double scalar)
{	
	// A = (A)+(s*B)
	/* code */
	int N = size;
	int elements = size*size;
	for (size_t i=0; i<elements; i++){
		A[i] = A[i] + scalar * B[i];
	}
}

struct Matrix {
	int size;
	double matrix[100];
	// row first form 
	// Each matrix can have max size of 10x10
};

void solve_bbd( int nMat, struct Matrix * matA, struct Matrix * matB, struct Matrix * matC, struct Matrix * matG){

	// Initialize Add WorkSpace
	int size = matA[0].size;
	int elements = size*size;
	double AddG[elements];
	for(int i =0; i<elements; i++) AddG[i] = 0;

	double AddB[elements];
	for(int i =0; i<elements; i++) AddB[i] = 0;

	#pragma omp prallel for
	for (int i=0; i<nMat; i++){

		printf("\n #%d ", i);
		printf("%s\n", "Matrix Before Inversion");
		print_matrix(matA[i].size, matA[i].matrix);

		// B = inv(A)*C*B
		// G = inv(A)*C*G
		matrix_inv_m(matA[i].size, matA[i].matrix, matC[i].matrix, matB[i].matrix, matG[i].matrix);
		// AddG = AddG + G 
		matrix_add(matA[i].size, AddG, matG[i].matrix, 1.0);
		// AddB = AddB + B 
		matrix_add(matA[i].size, AddB, matB[i].matrix, 1.0);

		printf("#%d ", i);
		printf("%s\n", "Matrix After Inversion");
		print_matrix(matG[i].size, matG[i].matrix);

	}

	printf("%s\n", "AddG Matrix WorkSpace After Addition");
	print_matrix(size, AddG);

	printf("%s\n", "AddB Matrix WorkSpace After Addition");
	print_matrix(size, AddB);

	printf("%s\n", "matA[nMat]");
	print_matrix(size, matA[nMat].matrix);

	printf("%s\n", "matG[nMat]");
	print_matrix(size, matG[nMat].matrix);

	matrix_add(size, matG[nMat].matrix, AddG, -1.0);
	matrix_add(size, matA[nMat].matrix, AddB, -1.0);

	printf("%s\n", "After sub matA[nMat]");
	print_matrix(size, matA[nMat].matrix);

	printf("%s\n", "After sub matG[nMat]");
	print_matrix(size, matG[nMat].matrix);

	matrix_invm(size, matA[nMat].matrix, matG[nMat].matrix, AddG);

	printf("%s\n", "After final Inversion matA[nMat]");
	print_matrix(size, AddG);

}

int main() {

	int nMat = 2;
	struct Matrix matA[nMat+1];
	struct Matrix matB[nMat];
	struct Matrix matC[nMat];
	struct Matrix matG[nMat+1];

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

	solve_bbd(nMat, matA, matB, matC, matG);

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

