

struct Matrix {
	// row first form 
	// Each matrix can have max size of 10x10
	int nrow;
	int ncol;
	double matval[10*10];
};

void solve_bbd( int nMat, int m, int n,
				struct Matrix * matA, struct Matrix * matB, 
				struct Matrix * matC, struct Matrix * matG, 
				struct Matrix * matX, struct Matrix  matAN, 
				struct Matrix  matGN, struct Matrix  matXN );

void solve_bbd_full( int N, int m, int n,
				struct Matrix * matA, struct Matrix * matB, 
				struct Matrix * matC, struct Matrix * matG, 
				struct Matrix  matAN, 
				struct Matrix  matGN);


int print_matrix(int nrow, int ncol, double * M);

