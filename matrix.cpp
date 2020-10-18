#include <iostream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "matrix.hpp"


Matrix::Matrix(int _height, int _width) : 
	width(_width), height(_height) {}

Matrix::Matrix(const Matrix& rhs) :
	width(rhs.width),
	height(rhs.height)
	{
		for (unsigned int i = 0; i < height; i++) {
			std::copy(rhs.m_data[i], rhs.m_data[i] + width, m_data[i]);
		}
	}

Matrix::~Matrix() {
	// Reset
	width = 0;
	height = 0;
}

Matrix& Matrix::operator=(const Matrix& rhs) {
    Matrix temp(rhs);
    swap(temp);
    return *this;
}

Matrix Matrix::operator*(const Matrix& rhs){
		
	if (width != rhs.height)
		throw std::invalid_argument( "ERROR: Matrix dimensions do not correspond! Aborting multiplication.\n" );
	
	Matrix result(height, rhs.width);
	float sum = 0;

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < rhs.width; j++) {  //j = cols
			sum = 0;
			for (unsigned int k = 0; k < width; ++k) {
				sum += m_data[i][k] * rhs.m_data[k][j];
			}
			result.m_data[i][j] = sum;
		}
	}
	return result;
}

Vertex Matrix::operator*(const Vertex& rhs) {

	Vertex result;
	float sum;

	for (unsigned int i = 0; i < height; i++) {
		sum = 0;
		for (unsigned int j = 0; j < width; j++) {
			 sum += m_data[i][j] * rhs.m_data[j];
		}
		result.m_data[i] = sum;
	}

	return result;
}

void Matrix::initData(const float *m, bool isVec) {

	int counter = 0;

	unsigned int for_j_ceiling = width;
	if (isVec) for_j_ceiling = 1;

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < for_j_ceiling; j++) {  //j = cols
			m_data[i][j] = m[counter];
			++counter;
		}
	}
}

void Matrix::makeIdentity() {

	for (unsigned int i = 0; i < height; i++) {  //i = rows
		for (unsigned int j = 0; j < width; j++) {  //j = cols
			if (i == j) m_data[i][j] = 1;
			else		m_data[i][j] = 0;
		}
	}
}

void Matrix::swap(Matrix& rhs) {
	std::swap(m_data, rhs.m_data);
	std::swap(width, rhs.width);
	std::swap(height, rhs.height);
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