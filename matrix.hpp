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
	Matrix() : width(4), height(4) { };

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
	/// Computes recurstion of current matrix
	///	WARNING: black magic happening here
	///		@return inverse of the current matrix
	Matrix inverse() {
		float A2323 = m_data[2][2] * m_data[3][3] - m_data[2][3] * m_data[3][2];
		float A1323 = m_data[2][1] * m_data[3][3] - m_data[2][3] * m_data[3][1];
		float A1223 = m_data[2][1] * m_data[3][2] - m_data[2][2] * m_data[3][1];
		float A0323 = m_data[2][0] * m_data[3][3] - m_data[2][3] * m_data[3][0];
		float A0223 = m_data[2][0] * m_data[3][2] - m_data[2][2] * m_data[3][0];
		float A0123 = m_data[2][0] * m_data[3][1] - m_data[2][1] * m_data[3][0];
		float A2313 = m_data[1][2] * m_data[3][3] - m_data[1][3] * m_data[3][2];
		float A1313 = m_data[1][1] * m_data[3][3] - m_data[1][3] * m_data[3][1];
		float A1213 = m_data[1][1] * m_data[3][2] - m_data[1][2] * m_data[3][1];
		float A2312 = m_data[1][2] * m_data[2][3] - m_data[1][3] * m_data[2][2];
		float A1312 = m_data[1][1] * m_data[2][3] - m_data[1][3] * m_data[2][1];
		float A1212 = m_data[1][1] * m_data[2][2] - m_data[1][2] * m_data[2][1];
		float A0313 = m_data[1][0] * m_data[3][3] - m_data[1][3] * m_data[3][0];
		float A0213 = m_data[1][0] * m_data[3][2] - m_data[1][2] * m_data[3][0];
		float A0312 = m_data[1][0] * m_data[2][3] - m_data[1][3] * m_data[2][0];
		float A0212 = m_data[1][0] * m_data[2][2] - m_data[1][2] * m_data[2][0];
		float A0113 = m_data[1][0] * m_data[3][1] - m_data[1][1] * m_data[3][0];
		float A0112 = m_data[1][0] * m_data[2][1] - m_data[1][1] * m_data[2][0];

		float det = m_data[0][0] * (m_data[1][1] * A2323 - m_data[1][2] * A1323 + m_data[1][3] * A1223)
			- m_data[0][1] * (m_data[1][0] * A2323 - m_data[1][2] * A0323 + m_data[1][3] * A0223)
			+ m_data[0][2] * (m_data[1][0] * A1323 - m_data[1][1] * A0323 + m_data[1][3] * A0123)
			- m_data[0][3] * (m_data[1][0] * A1223 - m_data[1][1] * A0223 + m_data[1][2] * A0123);
		det = 1 / det;

		Matrix ret;
		ret.m_data[0][0] = det * (m_data[1][1] * A2323 - m_data[1][2] * A1323 + m_data[1][3] * A1223);
		ret.m_data[0][1] = det * -(m_data[0][1] * A2323 - m_data[0][2] * A1323 + m_data[0][3] * A1223);
		ret.m_data[0][2] = det * (m_data[0][1] * A2313 - m_data[0][2] * A1313 + m_data[0][3] * A1213);
		ret.m_data[0][3] = det * -(m_data[0][1] * A2312 - m_data[0][2] * A1312 + m_data[0][3] * A1212);
		ret.m_data[1][0] = det * -(m_data[1][0] * A2323 - m_data[1][2] * A0323 + m_data[1][3] * A0223);
		ret.m_data[1][1] = det * (m_data[0][0] * A2323 - m_data[0][2] * A0323 + m_data[0][3] * A0223);
		ret.m_data[1][2] = det * -(m_data[0][0] * A2313 - m_data[0][2] * A0313 + m_data[0][3] * A0213);
		ret.m_data[1][3] = det * (m_data[0][0] * A2312 - m_data[0][2] * A0312 + m_data[0][3] * A0212);
		ret.m_data[2][0] = det * (m_data[1][0] * A1323 - m_data[1][1] * A0323 + m_data[1][3] * A0123);
		ret.m_data[2][1] = det * -(m_data[0][0] * A1323 - m_data[0][1] * A0323 + m_data[0][3] * A0123);
		ret.m_data[2][2] = det * (m_data[0][0] * A1313 - m_data[0][1] * A0313 + m_data[0][3] * A0113);
		ret.m_data[2][3] = det * -(m_data[0][0] * A1312 - m_data[0][1] * A0312 + m_data[0][3] * A0112);
		ret.m_data[3][0] = det * -(m_data[1][0] * A1223 - m_data[1][1] * A0223 + m_data[1][2] * A0123);
		ret.m_data[3][1] = det * (m_data[0][0] * A1223 - m_data[0][1] * A0223 + m_data[0][2] * A0123);
		ret.m_data[3][2] = det * -(m_data[0][0] * A1213 - m_data[0][1] * A0213 + m_data[0][2] * A0113);
		ret.m_data[3][3] = det * (m_data[0][0] * A1212 - m_data[0][1] * A0212 + m_data[0][2] * A0112);
		return ret;
	}

	/// Fills the matrix with given data
	///		@param m[in] matrix represented as an 16-float long 1D array
	void initData(const float *m) {

		int counter = 0;

		for (unsigned int i = 0; i < height; i++) {  //i = rows
			for (unsigned int j = 0; j < width; j++) {  //j = cols
				m_data[j][i] = m[counter];
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
	inline void swap(Matrix& rhs) {
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