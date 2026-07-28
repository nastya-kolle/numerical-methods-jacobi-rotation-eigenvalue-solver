#pragma once
#include "JacobiSolver.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include "LinearAlgebra.h"

struct ExperimentResult {
	int N;
	double lambda_min, lambda_max;
	double max_off_diag_final;
	int avg_iterations;
	double avg_lambda_error;
	double avg_accuracy_r;

	ExperimentResult() {
		N = 0;
		lambda_min = 0.0;
		lambda_max = 0.0;
		max_off_diag_final = 0.0;
		avg_iterations = 0;
		avg_lambda_error = 0.0;
		avg_accuracy_r = 0.0;
	}
};

class ComputationalExperiment {
public:

	// Запуск для ОДНОЙ пары (N, [λ_min, λ_max])
	static ExperimentResult run(int N, double lambda_min, double lambda_max);

	// Вспомогательные функции:
	static SquareMatrix generateOrthoronalMatrix(int N);
	static SquareMatrix generateSymmetricMatrix(int N, const Vector& eigen_values);

	// Вывод одного результата (формат таблицы — построчно)
	static void printHeader();
	static void printRow(int index, const ExperimentResult& res);
	static void printFooter();
};