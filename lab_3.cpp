#include "ComputationalExperiment.h"
#include "JacobiSolver.h"
#include "LinearAlgebra.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include <iostream>

using namespace std;
int main()
{
    setlocale(LC_ALL, "Russian");
    // Пример: sizes = {10, 30}, lambda_ranges = {(-2,2), (-50,50)}
    Vector sizes(2);
    sizes[1] = 10;
    sizes[2] = 30;

    Vector lambda_mins(2);
    Vector lambda_maxs(2);
    lambda_mins[1] = -2.0;
    lambda_maxs[1] = 2.0;
    lambda_mins[2] = -50.0;
    lambda_maxs[2] = 50.0;

    ComputationalExperiment::printHeader();

    int index = 1;
    for (int i = 1; i <= sizes.get_size(); i++) {
        int N = static_cast<int>(sizes[i]);
        for (int j = 1; j <= lambda_mins.get_size(); j++) {
            auto res = ComputationalExperiment::run(N, lambda_mins[j], lambda_maxs[j]);
            ComputationalExperiment::printRow(index++, res);
        }
    }
    ComputationalExperiment::printFooter();
}
