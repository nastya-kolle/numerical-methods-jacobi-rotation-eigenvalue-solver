#include "JacobiSolver.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;

JacobiResult JacobiSolver::solve(const SquareMatrix& A, double eps, int max_rotations)
{
	int N = A.get_size();
	if (N < 1) {
		throw invalid_argument("Matrix size must be at least 1");
	}

	SquareMatrix A0 = A;
	SquareMatrix V(N);
	V.fill_identity();

	int rotations = 0;

	while (rotations < max_rotations) {
		// --- Step 1: find largest off-diagonal element (p < q) ---
		int p = 0, q = 0;
		double max_elem = 0.0;
		for (int i = 1; i <= N; i++) {
			for (int j = i + 1; j <= N; j++) {
				double value = std::abs(A0[i][j]);
				if (value > max_elem) {
					max_elem = value;
					p = i;
					q = j;
				}
			}
		}
		// Convergence test
		if (max_elem < eps) {
			break;
		}
		// --- Step 2: compute rotation parameters ---
		double a_pp = A0[p][p];
		double a_qq = A0[q][q];
		double a_pq = A0[p][q];

		double cos_a, sin_a;
		if (std::abs(a_qq - a_pp) < 1e-15) {
			// Rotation by 45° (τ = 0 → θ = π/4)
			cos_a = sin_a = sqrt(0.5);
		}
		else {
			double ctg_2a = (a_qq - a_pp) / (2.0 * a_pq);
			double tg_a;
			if (ctg_2a >= 0) {
				tg_a = 1.0 / (ctg_2a + sqrt(ctg_2a * ctg_2a + 1.0));
			}
			else {
				tg_a = 1.0 / (ctg_2a - sqrt(ctg_2a * ctg_2a + 1.0));
			}
			cos_a = 1.0 / (sqrt(tg_a * tg_a + 1));
			sin_a = cos_a * tg_a;
		}

		// --- Step 3: update A0 ← J(p,q,θ)^T * A0 * J(p,q,θ) ---
		// Update only rows/cols p and q (symmetric)
		for (int k = 1; k <= N; k++) {
			if (k != p && k != q) {
				double a_pk = A0[p][k];
				double a_qk = A0[q][k];
				// A0[p][k] = A0[k][p] = c*apk - s*aqk;
				// A0[q][k] = A0[k][q] = s*apk + c*aqk;
				double new_a_pk= cos_a * a_pk - sin_a * a_qk;
				double new_a_qk = sin_a * a_pk + cos_a * a_qk;
				A0[p][k] = A0[k][p] = new_a_pk;
				A0[q][k] = A0[k][q] = new_a_qk;
			}
		}

		// Diagonal and (p,q) entry
		double new_a_pp = cos_a * cos_a * a_pp - 2.0 * cos_a * sin_a * a_pq + sin_a * sin_a * a_qq;
		double new_a_qq = sin_a * sin_a * a_pp + 2.0 * cos_a * sin_a * a_pq + cos_a * cos_a * a_qq;

		A0[p][p] = new_a_pp;
		A0[q][q] = new_a_qq;
		A0[p][q] = A0[q][p] = 0.0;

		// --- Step 4: update V ← V * J(p,q,θ) (right rotation) ---
		for (int i = 1; i <= N; i++) {
			double V_ip = V[i][p];
			double V_iq = V[i][q];
			double new_V_ip = cos_a * V_ip - sin_a * V_iq;
			double new_V_iq = sin_a * V_ip + cos_a * V_iq;
			V[i][p] = new_V_ip;
			V[i][q] = new_V_iq;
		}
		rotations++;
	}
	// --- Extract eigenvalues into Vector a (1-based!) ---
	Vector eigen_values(N);
	for (int i = 1; i <= N; i++) {
		eigen_values[i] = A0[i][i];
	}

	// --- Compute accuracy measure r = max_ij |(A·T − T·Λ)_ij| ---
	// Build Λ = diag(eigenvalues)
	SquareMatrix Lambda(N);
	Lambda.fill_zero();
	for (int i = 1; i <= N; i++) {
		Lambda[i][i] = eigen_values[i];
	}
	SquareMatrix AT = A * V; // A·T
	SquareMatrix TL = V * Lambda; // T·Λ
	SquareMatrix diff = AT - TL;

	double max_res = 0.0;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			double val = std::abs(diff[i][j]);
			if (val > max_res) max_res = val;
		}
	}

	double max_off = 0.0;
	for (int i = 1; i < N; i++) {
		for (int j = i + 1; j <= N; j++) {
			max_off = std::max(max_off, std::abs(A0[i][j]));
		}
	}

	// --- Assemble result ---
	JacobiResult res;
	res.IER = (rotations < max_rotations) ? 0 : rotations;
	res.T = std::move(V);
	res.a = std::move(eigen_values);
	res.K = rotations;
	res.r = max_res;
	res.final_max_off_diag = max_off;

	return res;
}
