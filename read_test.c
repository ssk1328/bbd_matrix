#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

// gcc t2.c -llapack -std=c99

int print_matrix(int nrow, int ncol, double * M)
{

	for (int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			printf("%f ", M[i*ncol + j]);
		}
		printf("\n");
	}
	return 0;
}


struct Matrix {
	// row first form 
	// Each matrix can have max size of 10x10
	// TODO Add support for having non square matrices
	int nrow;
	int ncol;
	double matval[10*10];
};

void solve_bbd( int nMat, struct Matrix * matA, struct Matrix * matB, struct Matrix * matC, struct 	Matrix * matG, struct 	Matrix * matX);

int main() {

// ---------------------------------------------------------------------
// END OF READ MATRIX DATA FROM FILE
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

	for (int j=0; j<N-1; j++){

		// Initialize A
		for	(int i=0; i<m*m; i++){
			fscanf(fp, "%lf", &matA[j].matval[i]);
		}
		print_matrix(matA[j].nrow, matA[j].ncol, matA[j].matval);

		// Initialize B
		for	(int i=0; i<m*n; i++){
			fscanf(fp, "%lf", &matB[j].matval[i]);
		}
		print_matrix(matB[j].nrow, matB[j].ncol, matB[j].matval);

		// Initialize C
		for	(int i=0; i<n*m; i++){
			fscanf(fp, "%lf", &matC[j].matval[i]);
		}
		print_matrix(matC[j].nrow, matC[j].ncol, matC[j].matval);

		// Initialize G
		for	(int i=0; i<m*1; i++){
			fscanf(fp, "%lf", &matG[j].matval[i]);
		}
		print_matrix(matG[j].nrow, matG[j].ncol, matG[j].matval);
	}

	// Initialize AN
	for	(int i=0; i<n*n; i++){
		fscanf(fp, "%lf", &matAN.matval[i]);
	}
	print_matrix(matAN.nrow, matAN.ncol, matAN.matval);

	// Initialize GN
	for	(int i=0; i<n*1; i++){
		fscanf(fp, "%lf", &matGN.matval[i]);
	}
	print_matrix(matGN.nrow, matGN.ncol, matGN.matval);

	fclose(fp);

// ---------------------------------------------------------------------
// END OF READ MATRIX DATA FROM FILE
// ---------------------------------------------------------------------

	printf("******************** Matrix Read Done ******************\n" );

	solve_bbd(N-1, matA, matB, matC, matG, matX, matAN, matGN, matXN);

	return 0;
}

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