#ifndef MATRIX_H
#define MATRIX_H

#include "vertex.h"

/// Class representing a 4x4 flat matrix, can be used as a representation of points
class Matrix {
public:
	/* --- VARIABLES --- */
	/// Width of the matrix
	unsigned int width;
	/// Height of the matrix
	unsigned int height;
	/// The matrix itself
	float m_data[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

	/* --- CONSTRUCTORS & destructors ---*/
	/// Custom Matrix constructor
	Matrix(int _height, int _width) :
		width(_width), height(_height) {}

	/// Default matrix constructor - we are using 4x4 matrices
	Matrix() : width(4), height(4) {};

	/// Matrix Copy constructor
	Matrix(const Matrix& rhs) :
		width(rhs.width),
		height(rhs.height)
	{
		for (unsigned int i = 0; i < height; i++) {
			std::copy(rhs.m_data[i], rhs.m_data[i] + width, m_data[i]);
		}
	}

	/// Matrix Destructor
	~Matrix() {
		width = 0;
		height = 0;
	}

	/* --- OPERATORS ---*/
	/// Operator * for matrix mulitplication with another matrix
	///		@return resulting matrix
	Matrix operator*(const Matrix& rhs) {

		/*if (width != rhs.height)
			throw std::invalid_argument("ERROR: Matrix dimensions do not correspond! Aborting multiplication.\n");*/

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

	/// Operator * for matrix mulitplication with a vector
	///		@param rhs[in] which vertex should be multiplied with this matrix
	///		@return resulting vector
	Vertex operator*(const Vertex& rhs) {

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

	/// Matrix assignment Operator
	///		@param rhs[in] which value to assign
	Matrix& operator=(const Matrix& rhs) {
		Matrix temp(rhs);
		swap(temp);
		return *this;
	}

	/* --- FUNCTIONS ---*/
	/// Fills the matrix with given data
	///		@param m[in] matrix represented as an 16-float long 1D array
	void initData(const float *m) {

		int counter = 0;

		for (unsigned int i = 0; i < height; i++) {  //i = rows
			for (unsigned int j = 0; j < width; j++) {  //j = cols
				m_data[i][j] = m[counter];
				++counter;
			}
		}
	}

	/// Makes the current matrix an identity matrix
	void makeIdentity() {

		for (unsigned int i = 0; i < height; i++) {  //i = rows
			for (unsigned int j = 0; j < width; j++) {  //j = cols
				if (i == j) m_data[i][j] = 1;
				else		m_data[i][j] = 0;
			}
		}
	}

	/// Function for move constructor
	///		@param rhs[in] reference to the matrix with which data should be swapped 
	void swap(Matrix& rhs) {
		std::swap(m_data, rhs.m_data);
		std::swap(width, rhs.width);
		std::swap(height, rhs.height);
	}

	/// Outputs the matrix to the console
	void print() {
		for (unsigned int i = 0; i < height; ++i) {
			std::cout << "[";
			for (unsigned int j = 0; j < width; ++j) {
				std::cout << " " << m_data[i][j];
			}
			std::cout << " ]\n";
		}
	}
};

#endif