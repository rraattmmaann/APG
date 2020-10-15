#pragma once

#include <iostream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <thread>

#include "matrix.hpp"

// Time to miliseconds, for messuring commputation time
template <typename TimePoint>
std::chrono::milliseconds to_ms(TimePoint tp) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(tp);
}

// Custom constructor
Matrix::Matrix(int _height, int _width) : 
	width(_width), height(_height) {

	allocate();

	/*for (unsigned int i = 0; i < height; ++i) {
		for (unsigned int j = 0; j < width; ++j) {
			std::cin >> m_data[i][j];
		}
	}*/
}

// Copy constructor
Matrix::Matrix(const Matrix& rhs) :
	width(rhs.width),
	height(rhs.height)
	{
		allocate();
		for (unsigned int i = 0; i < height; i++) {
			std::copy(rhs.m_data[i], rhs.m_data[i] + width, m_data[i]);
		}
	}

// Destructor
Matrix::~Matrix() {
	if (m_data == NULL) return;

	// Free the memory.
	for (unsigned int i = 0; i < height; ++i ) {
		delete[] m_data[i];
	}
	delete[] m_data;

	// Reset
	width = 0;
	height = 0;
	m_data = NULL;
}

// Assignment Operator
Matrix& Matrix::operator=(const Matrix& rhs) {
    Matrix temp(rhs);
    swap(temp);
    return *this;
}

// Allcoate matrix
void Matrix::allocate() {

	m_data = new float*[height];
	for (unsigned int i = 0; i < height; i++ ) {
		m_data[i] = new float[width]();
	}
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
	Matrix b = rhs.transpose();

	result = multiplySequential(*this, b, result);

	return result;
}

// Transpose this matrix and return transposed matrix
Matrix Matrix::transpose() const {

    Matrix transposed(width, height);

    for (unsigned int i = 0; i < height; i++) {
        for (unsigned int j = 0; j < width; j++) {
            transposed.m_data[j][i] = m_data[i][j];
        }
    }

    return transposed;
}

// Sequential multiplication of two matrices
Matrix& Matrix::multiplySequential(const Matrix& a, const Matrix& b, Matrix& result) {

	float sum = 0;
	//auto start = std::chrono::high_resolution_clock::now();

	for (unsigned int i=0; i < a.height; i++) {  //i = rows
		for (unsigned int j=0; j < b.height; j++) {  //j = cols
			sum = 0;
			for(unsigned int k = 0; k < a.width; ++k) {
                sum += a.m_data[i][k] * b.m_data[j][k];
            }
            result.m_data[i][j] = sum;
		}
	}
	//auto end = std::chrono::high_resolution_clock::now();

	//std::cout << "Sequential multiplication took " << to_ms(end - start).count() << "ms.\n";

	return result;
}

void Matrix::makeIdentity() {

	if (m_data == NULL) allocate();

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < width; j++) {  //j = cols
			if (i == j) m_data[i][j] = 1;
			else		m_data[i][j] = 0;
		}
	}
}

void Matrix::initData(const float *m) {

	if (m_data == NULL) allocate();

	int counter = 0;
	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < width; j++) {  //j = cols
			m_data[j][i] = m[counter];
			++counter;
		}
	}
}