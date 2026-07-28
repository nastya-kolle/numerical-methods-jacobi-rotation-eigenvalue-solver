#pragma once
#include "SquareMatrix.h"
#include "Vector.h"

struct JacobiResult {
	int IER; // 0 — success; >0 — max rotations reached (K = IER)
	SquareMatrix T; // N×N matrix: column i is eigenvector for a[i]
	Vector a; // eigenvalues: a[i], i = 1..N
	int K; // number of rotations performed
	double r; // accuracy: max |A·T − T·Λ|_ij
	double final_max_off_diag; // max |A0[i][j]| for i≠j after iterations
};

class JacobiSolver {
public:
	static JacobiResult solve(const SquareMatrix& A, double eps = 1e-8, int max_rotations = 1000);
};
