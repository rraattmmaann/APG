#include <iostream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "matrix.hpp"

// Custom constructor
Matrix::Matrix(int _height, int _width) : 
	width(_width), height(_height) {}

// Copy constructor
Matrix::Matrix(const Matrix& rhs) :
	width(rhs.width),
	height(rhs.height)
	{
		for (unsigned int i = 0; i < height; i++) {
			std::copy(rhs.m_data[i], rhs.m_data[i] + width, m_data[i]);
		}
	}

// Destructor
Matrix::~Matrix() {
	if (m_data == NULL) return;

	// Reset
	width = 0;
	height = 0;
}

// Assignment Operator
Matrix& Matrix::operator=(const Matrix& rhs) {
    Matrix temp(rhs);
    swap(temp);
    return *this;
}

// Function for move constructor
void Matrix::swap(Matrix& rhs) {
    std::swap(m_data, rhs.m_data);
    std::swap(width, rhs.width);
    std::swap(height, rhs.height);
}

// Operator * for matrix mulitplication
Matrix Matrix::operator*(const Matrix& rhs){
		
	if (width != rhs.height)
		throw std::invalid_argument( "ERROR: Matrix dimensions do not correspond! Aborting multiplication.\n" );
	
	Matrix result(height, rhs.width);

	return multiplySequential(*this, rhs, result);
}

Matrix Matrix::operator*(const float& rhs) {

	Matrix result(height, width);
	
	// In this case we only use this operator to multiply a vector by a value,
	// so we do not need to multiply the whole matrix
	for (unsigned int i = 0; i < height; i++) {
		for (unsigned int j = 0; j < 1; j++) {				
			result.m_data[i][j] = m_data[i][j] *= rhs;
		}
	}

	return result;
}

// Sequential multiplication of two matrices
Matrix& Matrix::multiplySequential(const Matrix& a, const Matrix& b, Matrix& result) {

	float sum = 0;

	for (unsigned int i=0; i < a.height; i++) {  //i = rows
		for (unsigned int j=0; j < b.width; j++) {  //j = cols
			sum = 0;
			for(unsigned int k = 0; k < a.width; ++k) {
                sum += a.m_data[i][k] * b.m_data[k][j];
            }
            result.m_data[i][j] = sum;
		}
	}

	return result;
}

void Matrix::makeIdentity() {

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < width; j++) {  //j = cols
			if (i == j) m_data[i][j] = 1;
			else		m_data[i][j] = 0;
		}
	}
}

void Matrix::initData(const float *m, bool isVec) {

	int counter = 0;

	int for_j_ceiling = width;
	if (isVec) for_j_ceiling = 1;

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < for_j_ceiling; j++) {  //j = cols
			m_data[i][j] = m[counter];
			++counter;
		}
	}
}

void Matrix::print() {
	for (unsigned int i = 0; i < height; ++i) {
		std::cout << "[";
		for (unsigned int j = 0; j < width; ++j) {
			std::cout << " " << m_data[i][j];
		}
		std::cout << " ]\n";
	}
}