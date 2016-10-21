#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include "bbdf.h"


// gcc t2.c -llapack -std=c99

int main() {

    const rlim_t kStackSize = 1024 * 1024 * 1024;   // min stack size = 32 MB
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

// ---------------------------------------------------------------------
// START OF READ MATRIX DATA FROM FILE
// ---------------------------------------------------------------------



	FILE *fp;
	fp = fopen("data.txt", "r");

	int N;
	int m, n;

	fscanf(fp, "%d %d %d", &N, &m, &n );

	struct Matrix matA[N-1];
	struct Matrix matB[N-1];
	struct Matrix matC[N-1];
	struct Matrix matG[N-1];
	struct Matrix matX[N-1];
	struct Matrix matAN;
	struct Matrix matGN;
	struct Matrix matXN;

	matAN.nrow = n;
	matAN.ncol = n;
	matGN.nrow = n;
	matGN.ncol = 1;
	matXN.nrow = n;
	matXN.ncol = 1;

	for (int i=0; i<N-1; i++){
		matA[i].nrow = m;
		matA[i].ncol = m;
		matB[i].nrow = m;
		matB[i].ncol = n;
		matC[i].nrow = n;
		matC[i].ncol = m;
		matG[i].nrow = m;
		matG[i].ncol = 1;
		matX[i].nrow = m;
		matX[i].ncol = 1;
	}

	for (int j=0; j<N-1; j++){

		// Initialize A
		for	(int i=0; i<m*m; i++){
			fscanf(fp, "%lf", &matA[j].matval[i]);
		}
//		print_matrix(matA[j].nrow, matA[j].ncol, matA[j].matval);

		// Initialize B
		for	(int i=0; i<m*n; i++){
			fscanf(fp, "%lf", &matB[j].matval[i]);
		}
//		print_matrix(matB[j].nrow, matB[j].ncol, matB[j].matval);

		// Initialize C
		for	(int i=0; i<n*m; i++){
			fscanf(fp, "%lf", &matC[j].matval[i]);
		}
//		print_matrix(matC[j].nrow, matC[j].ncol, matC[j].matval);

		// Initialize G
		for	(int i=0; i<m*1; i++){
			fscanf(fp, "%lf", &matG[j].matval[i]);
		}
//		print_matrix(matG[j].nrow, matG[j].ncol, matG[j].matval);
	}

	// Initialize AN
	for	(int i=0; i<n*n; i++){
		fscanf(fp, "%lf", &matAN.matval[i]);
	}
//	print_matrix(matAN.nrow, matAN.ncol, matAN.matval);

	// Initialize GN
	for	(int i=0; i<n*1; i++){
		fscanf(fp, "%lf", &matGN.matval[i]);
	}
//	print_matrix(matGN.nrow, matGN.ncol, matGN.matval);

	fclose(fp);

// ---------------------------------------------------------------------
// END OF READ MATRIX DATA FROM FILE
// ---------------------------------------------------------------------

	printf("******************** Matrix Read Done ******************\n" );

	clock_t tStart = clock();
//	solve_bbd( N-1 , m, n, matA, matB, matC, matG, matX, matAN, matGN, matXN);
	clock_t end = clock();
	printf("Total Time taken: %.4fms\n", (double)(end - tStart)/(CLOCKS_PER_SEC*0.001));

	solve_bbd_full(N , m, n, matA, matB, matC, matG, matAN, matGN);

	return 0;
}
