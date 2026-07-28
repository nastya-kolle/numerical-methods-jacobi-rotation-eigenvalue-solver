#include <string>
#include <utility>
#include "Vector.h"
#include "SquareMatrix.h"
#include "LinearAlgebra.h"

using namespace std;

// Ортогонализация вектора v относительно u
void LinearAlgebra::orthogonalizeAgainst(Vector& v, const Vector& u)
{
	if (u.get_size() != v.get_size())
		throw std::invalid_argument("Vectors must have same length");

	double dot = 0.0;
	for (int i = 1; i <= u.get_size(); ++i)
		dot += u[i] * v[i];

	for (int i = 1; i <= v.get_size(); ++i)
		v[i] -= dot * u[i];
}

// Нормировка вектора
bool LinearAlgebra::normalize(Vector& v, double epsilon )
{
	double norm_sq = 0.0;
	for (int i = 1; i <= v.get_size(); ++i)
		norm_sq += v[i] * v[i];

	if (norm_sq < epsilon * epsilon)
		return false;

	double norm = sqrt(norm_sq);
	for (int i = 1; i <= v.get_size(); ++i)
		v[i] /= norm;

	return true;
}

