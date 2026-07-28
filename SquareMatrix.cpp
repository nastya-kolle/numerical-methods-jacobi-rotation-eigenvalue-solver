#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include "Vector.h"
#include "SquareMatrix.h"

using namespace std;

SquareMatrix::SquareMatrix(int size_) : size(size_)
{
	if (size_ < 1) throw invalid_argument("Matrix size must be at least 1");
	rows = new Vector[size];
	for (int i = 0; i < size_; i++) {
		rows[i] = Vector(size);
	}
}
SquareMatrix::SquareMatrix(const SquareMatrix& other) : size(other.size)
{
	rows = new Vector[size];
	for (int i = 0; i < size; i++) {
		rows[i] = other.rows[i];
	}
}

SquareMatrix::~SquareMatrix()
{
	delete[] rows;
}

int SquareMatrix::get_size() const
{
	return size;
}

SquareMatrix& SquareMatrix::operator=(const SquareMatrix& other)
{
	if (this != &other) {
		delete[] rows;
		size = other.size;
		rows = new Vector[size];
		for (int i = 0; i < size; i++) {
			rows[i] = other.rows[i];
		}
	}
	return *this;
}

Vector& SquareMatrix::operator[](int index)
{
	if (index < 1 || index > size) {
		throw out_of_range("Row index out of range" + to_string(index));
	}
	return rows[index - 1];
}

const Vector& SquareMatrix::operator[](int index) const
{
	if (index < 1 || index > size) {
		throw out_of_range("Row index out of range" + to_string(index));
	}
	return rows[index - 1];
}

double& SquareMatrix::at(int i, int j)
{
	validate_indices(i, j);
	return rows[i - 1][j];
}

const double& SquareMatrix::at(int i, int j) const
{
	validate_indices(i, j);
	return rows[i - 1][j];
}

// Умножение на вектор: y = A * x
Vector SquareMatrix::operator*(const Vector& vec) const
{
	if (size != vec.get_size()) {
		throw invalid_argument("Matrix and vector sizes mismatch");
	}
	Vector res(size);

#pragma omp parallel for
	for (int i = 1; i <= size; i++) {
		double sum = 0.0;
		for (int j = 1; j <= size; j++) {
			sum += (*this)[i][j] * vec[j];
		}
		res[i] = sum;
	}
	return res;
}

SquareMatrix SquareMatrix::operator-(const SquareMatrix& matrix) const
{
	if (size != matrix.size) {
		throw invalid_argument("Matrix sizes must match for subraction");
	}
	SquareMatrix res(size);
	for (int i = 1; i <= size; i++) {
		res[i] = (*this)[i] - matrix[i];
	}
	return res;
}

SquareMatrix SquareMatrix::operator+(const SquareMatrix& matrix) const
{
	if (size != matrix.size) {
		throw invalid_argument("Matrix sizes must match for addition");
	}
	SquareMatrix res(size);
	for (int i = 1; i <= size; i++) {
		res[i] = (*this)[i] + matrix[i];
	}
	return res;
}

// Умножение матриц: C = A * B
SquareMatrix SquareMatrix::operator*(const SquareMatrix& matrix) const
{
	if (size != matrix.size) {
		throw invalid_argument("Matrix sizes must match for multiplication");
	}
	SquareMatrix res(size);
	for (int i = 1; i <= size; i++) {
		for (int j = 1; j <= size; j++) {
			double sum = 0.0;
			for (int k = 1; k <= size; k++) {
				sum += (*this)[i][k] * matrix[k][j];
			}
			res[i][j] = sum;
		}
	}
	return res;
}

SquareMatrix& SquareMatrix::operator-=(const SquareMatrix& other)
{
	if (size != other.size) {
		throw invalid_argument("Matrix sizes must match for subraction");
	}
	for (int i = 0; i < size; i++) {
		(*this)[i + 1] -= other[i + 1];
	}
	return *this;
}

SquareMatrix& SquareMatrix::operator+=(const SquareMatrix& other)
{
	if (size != other.size) {
		throw invalid_argument("Matrix sizes must match for addition");
	}
	for (int i = 0; i < size; i++) {
		(*this)[i + 1] += other[i + 1];
	}
	return *this;
}

SquareMatrix SquareMatrix::operator*(double scalar) const
{
	SquareMatrix res(size);
	for (int i = 0; i < size; i++) {
		res[i] = (*this)[i + 1] * scalar;
	}
	return res;
}

SquareMatrix& SquareMatrix::operator*=(double scalar)
{
	for (int i = 0; i < size; i++) {
		(*this)[i + 1] *= scalar;
	}
	return *this;
}

bool SquareMatrix::is_symmetric(double eps) const
{
	for (int i = 1; i <= size; i++) {
		for (int j = 1; j <= size; j++) {
			if (abs((*this)[i][j]) - (*this)[j][i] > eps)
				return false;
		}
	}
	return true;
}

SquareMatrix SquareMatrix::transpose() const
{
	SquareMatrix res(size);
	for (int i = 1; i <= size; i++) {
		for (int j = 1; j <= size; j++) {
			res[i][j] = (*this)[j][i];
		}
	}
	return res;
}

void SquareMatrix::fill_random(double min_val, double max_val)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(min_val, max_val);

	for (int i = 1; i <= size; i++)
		for (int j = 1; j <= size; j++)
			(*this)[i][j] = dist(gen);
}

void SquareMatrix::fill_identity()
{
	fill_zero();
	for (int i = 1; i <= size; i++)
		(*this)[i][i] = 1.0;
}

void SquareMatrix::fill_zero()
{
	for (int i = 1; i <= size; i++)
		for (int j = 1; j <= size; j++)
			(*this)[i][j] = 0.0;
}

void SquareMatrix::print() const
{
	cout << toString();
}

std::string SquareMatrix::toString() const
{
	std::ostringstream oss;
	for (int i = 1; i <= size; i++) {
		oss << "[";
		for (int j = 1; j <= size; j++) {
			oss << std::setw(10) << std::setprecision(4) << (*this)[i][j];
			if (j < size) oss << ", ";
		}
		oss << "]" << endl;
	}
	return oss.str();
}

void SquareMatrix::print_to_file(std::ofstream& file) const
{
	file << "Matrix " << size << "x" << size << ":\n";
	file << toString();
}

SquareMatrix operator*(double scalar, const SquareMatrix& matrix)
{
	return matrix * scalar;
}
