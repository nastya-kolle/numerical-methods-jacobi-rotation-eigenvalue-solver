#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <string>
#include "Vector.h"

using namespace std;
class SquareMatrix {
private:
	int size;
	Vector* rows;

	void validate_indices(int i, int j) const {
		if (i < 1 || i > size || j < 1 || j > size) {
			throw out_of_range("Matrix indices out of range (" + to_string(i) + ", " + to_string(j) + ")");
		}
	}

public:
	SquareMatrix(int size_ = 2);
	SquareMatrix(const SquareMatrix& other);
	~SquareMatrix();

	int get_size() const;

	// Присваивание
	SquareMatrix& operator=(const SquareMatrix& other);

	Vector& operator[](int index);
	const Vector& operator[](int index) const;

	double& at(int i, int j);
	const double& at(int i, int j) const;

	// Арифметика
	Vector operator*(const Vector& vec) const;
	SquareMatrix operator-(const SquareMatrix& matrix) const;
	SquareMatrix operator+(const SquareMatrix& matrix) const;
	SquareMatrix operator*(const SquareMatrix& matrix) const;
	SquareMatrix& operator-=(const SquareMatrix& other);
	SquareMatrix& operator+=(const SquareMatrix& other);

	// Умножение на скаляр
	SquareMatrix operator*(double scalar) const;
	SquareMatrix& operator*=(double scalar);
	friend SquareMatrix operator*(double scalar, const SquareMatrix& matrix);


	bool is_symmetric(double eps = 1e-9) const;
	SquareMatrix transpose() const;

	// Заполнение
	void fill_random(double min_val, double max_val);
	void fill_identity();
	void fill_zero();

	// Ввод-вывод
	void print() const;
	std::string toString() const;
	void print_to_file(std::ofstream& file) const;
};


















