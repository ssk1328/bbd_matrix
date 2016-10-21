#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>

// gcc t2.c -llapack -std=c99

// lapack and blas functions used

// Two step process to compute inverse
extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV, int * INFO);
extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO);

int print_matrix(int nrow, int ncol, double * M, FILE * fp)
{
	for (int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			fprintf(fp, "%g ", M[i*ncol + j]);
		}
		fprintf(fp, "\n");
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
	int nrow;
	int ncol;
	double matval[10*10];
};

int main() {

    const rlim_t kStackSize = 1024 * 1024 * 1024;   // min stack size = 1024 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }

	int N = 10000; 	// number of blocks in diagonal
	int m = 4;	// mxm is the size of Ai
	int n = 4;	// nxn is the size of AN
	FILE *fp;
	
	fp = fopen("data.txt", "w+");

	fprintf(fp, "%d\n", N);
	fprintf(fp, "%d %d\n", m, n);

	time_t t;
	srand((unsigned) time(&t));
	
	struct Matrix matA[N-1];
	struct Matrix matB[N-1];
	struct Matrix matC[N-1];
	struct Matrix matG[N-1];
	struct Matrix matAN;
	struct Matrix matGN;

	matAN.nrow = n;
	matAN.ncol = n;
	matGN.nrow = n;
	matGN.ncol = 1;

	for (int i=0; i<N-1; i++){
		matA[i].nrow = m;
		matA[i].ncol = m;
		matB[i].nrow = m;
		matB[i].ncol = n;
		matC[i].nrow = n;
		matC[i].ncol = m;
		matG[i].nrow = m;
		matG[i].ncol = 1;
	}

	int val_max = 1024;

	for (int j=0; j<N-1; j++){
		// Initialize A
		for	(int i=0; i<m*m; i++){
			matA[j].matval[i] = (rand() % val_max) ;
		}
		matrix_inv_m (matA[j].nrow, matA[j].matval) ;
		for	(int i=0; i<m*m; i++){
			matA[j].matval[i] = matA[j].matval[i]*1024*128;
		}
//		fprintf(fp, "A%d\n", j);
		print_matrix(matA[j].nrow, matA[j].ncol, matA[j].matval, fp);		

		// Initialize B
		for	(int i=0; i<m*n; i++){
			matB[j].matval[i] = rand() % val_max ;
		}
//		fprintf(fp, "B%d\n", j);
		print_matrix(matB[j].nrow, matB[j].ncol, matB[j].matval, fp);		

		// Initialize C
		for	(int i=0; i<n*m; i++){
			matC[j].matval[i] = rand() % val_max ;
		}
//		fprintf(fp, "C%d\n", j);
		print_matrix(matC[j].nrow, matC[j].ncol, matC[j].matval, fp);		
		
		// Initialize G
		for	(int i=0; i<m*1; i++){
			matG[j].matval[i] = rand() % val_max ;		
		}
//		fprintf(fp, "G%d\n", j);
		print_matrix(matG[j].nrow, matG[j].ncol, matG[j].matval, fp);		
	
	}

	// Initialize AN
	for	(int i=0; i<n*n; i++){
		matAN.matval[i] = rand() % val_max ;
	}
//	fprintf(fp, "AN\n");
	print_matrix(matAN.nrow, matAN.ncol, matAN.matval, fp);		

	// Initialize GN
	for	(int i=0; i<n*1; i++){
		matGN.matval[i] = rand() % val_max ;
	}
//	fprintf(fp, "GN\n");
	print_matrix(matGN.nrow, matGN.ncol, matGN.matval, fp);

	fclose(fp);
    printf("Write Done\n");

	return 0;
}

