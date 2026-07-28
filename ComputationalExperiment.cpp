#include "ComputationalExperiment.h"
#include "JacobiSolver.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include "LinearAlgebra.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>

ExperimentResult ComputationalExperiment::run(int N, double lambda_min, double lambda_max)
{
	ExperimentResult res;
	res.N = N;;
	res.lambda_min = lambda_min;
	res.lambda_max = lambda_max;

	std::vector<double > epsilons = { 1e-5, 1e-7, 1e-9 };
	int total_rotations = 0;
	double total_r = 0.0;
	double max_lambda_error = 0.0;

	for (double eps : epsilons) {
		Vector true_vals(N);
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dist(lambda_min, lambda_max);

		for (int i = 1; i <= N; i++) {
			true_vals[i] = dist(gen);
		}

		SquareMatrix A = generateSymmetricMatrix(N, true_vals);
		//cout << endl;
		auto jacobi_result = JacobiSolver::solve(A, eps, 10000);

		total_rotations += jacobi_result.K;
		total_r += jacobi_result.r;

		// Максимальная ошибка по λ
		vector<double> true_vec(N);
		vector<double> computed_vec(N);

		for (int i = 1; i <= N; i++) {
			true_vec[i - 1] = true_vals[i];
			computed_vec[i - 1] = jacobi_result.a[i];
		}

		std::sort(true_vec.begin(), true_vec.end());
		std::sort(computed_vec.begin(), computed_vec.end());

		double max_err = 0.0;
		for (int i = 0; i < N; i++) {
			double err = std::abs(true_vec[i] - computed_vec[i]);
			if (err > max_err) max_err = err;
		}
		if (max_err > max_lambda_error) max_lambda_error = max_err;

		if (eps < res.max_off_diag_final || res.max_off_diag_final == 0.0)
			res.max_off_diag_final = eps;
	}
	res.avg_iterations = total_rotations / epsilons.size();
	res.avg_lambda_error = max_lambda_error;
	res.avg_accuracy_r = total_r / epsilons.size();

	return res;
}

SquareMatrix ComputationalExperiment::generateOrthoronalMatrix(int N)
{
	SquareMatrix Q(N);
	LinearAlgebra la;
	Q.fill_random(-1.0, 1.0);

	for (int j = 1; j <= N; j++) {
		Vector col = Q[j];

		// --- Первый проход: модифицированный Грам–Шмидт ---
		for (int k = 1; k < j; k++) {
			// Проекция col на Q[k]
			double proj = 0.0;
			for (int i = 1; i <= N; i++) {
				proj += col[i] * Q[i][k];  // col · Q[:,k]
			}
			// Вычитаем проекцию: col ← col - proj * Q[:,k]
			for (int i = 1; i <= N; i++) {
				col[i] -= proj * Q[i][k];
			}
		}

		// --- Повторная ортогонализация (второй проход) ---
		for (int k = 1; k < j; k++) {
			double proj = 0.0;
			for (int i = 1; i <= N; i++) {
				proj += col[i] * Q[i][k];
			}
			// Повторно вычитаем — даже если proj мала
			for (int i = 1; i <= N; i++) {
				col[i] -= proj * Q[i][k];
			}
		}

		// Нормировка
		if (!la.normalize(col, 1e-12)) {
			throw std::runtime_error("Cannot normalize column " + std::to_string(j));
		}

		// Записываем обновлённый столбец
		for (int i = 1; i <= N; i++) {
			Q[i][j] = col[i];
		}
	}
	return Q;
}

SquareMatrix ComputationalExperiment::generateSymmetricMatrix(int N, const Vector& eigen_values)
{
	SquareMatrix V = generateOrthoronalMatrix(N);
	SquareMatrix Lambda(N);
	Lambda.fill_zero();
	for (int i = 1; i <= N; i++) {
		Lambda[i][i] = eigen_values[i];
	}

	return V * Lambda * V.transpose();
}

void ComputationalExperiment::printHeader()
{
	std::cout << "\n" << std::setw(80) << std::setfill('=') << "" << std::setfill(' ') << "\n";
	std::cout << "               РЕЗУЛЬТАТЫ ВЫЧИСЛИТЕЛЬНЫХ ЭКСПЕРИМЕНТОВ\n";
	std::cout << "               (метод вращений Якоби для симметричных матриц)\n";
	std::cout << std::setw(80) << std::setfill('=') << "" << std::setfill(' ') << "\n";

	std::cout << std::left
		<< std::setw(4) << "№"
		<< std::setw(8) << "N"
		<< std::setw(15) << "Диапазон l "
		<< std::setw(18) << "Макс. |A_ij|"
		<< std::setw(12) << "Ср. итер."
		<< std::setw(15) << "Ср. оценка l"
		<< std::setw(15) << "Ср. мера r"
		<< "\n";

	std::cout << std::setw(80) << std::setfill('-') << "" << std::setfill(' ') << "\n";
}

void ComputationalExperiment::printRow(int index, const ExperimentResult& res)
{
	std::ostringstream range;
	range << "[" << res.lambda_min << ";" << res.lambda_max << "]";

	std::cout << std::left
		<< std::setw(4) << index
		<< std::setw(8) << res.N
		<< std::setw(15) << range.str()
		<< std::setw(18) << std::scientific << std::setprecision(2) << res.max_off_diag_final
		<< std::setw(12) << res.avg_iterations
		<< std::setw(15) << std::fixed << std::setprecision(6) << res.avg_lambda_error
		<< std::setw(15) << std::scientific << std::setprecision(2) << res.avg_accuracy_r
		<< "\n";
}

void ComputationalExperiment::printFooter()
{
	std::cout << std::setw(80) << std::setfill('=') << "" << std::setfill(' ') << "\n";
}
