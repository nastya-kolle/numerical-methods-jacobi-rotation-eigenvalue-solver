#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include "Vector.h"
#include "SquareMatrix.h"

using namespace std;

class LinearAlgebra {

public:
	void orthogonalizeAgainst(Vector& v, const Vector& u);
	bool normalize(Vector& v, double epsilon = 1e-15);
	double rayleighQuotient(const SquareMatrix& A, const Vector& x);
	double computeAccuracyMeasure(const SquareMatrix& A, double lambda, const Vector& x);
	double cosineAngle(const Vector& v1, const Vector& v2);

};