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
	Matrix(int _height, int _width);

	/// Default matrix constructor - we are using 4x4 matrices
	Matrix() : width(4), height(4) {};

	/// Matrix Copy constructor
	Matrix(const Matrix& rhs);	

	/// Matrix Destructor
	~Matrix();

	/* --- OPERATORS ---*/
	/// Operator * for matrix mulitplication with another matrix
	///		@return resulting matrix
	Matrix operator*(const Matrix& rhs);

	/// Operator * for matrix mulitplication with a vector
	///		@param rhs[in] which vertex should be multiplied with this matrix
	///		@return resulting vector
	Vertex operator*(const Vertex& rhs);

	/// Matrix assignment Operator
	///		@param rhs[in] which value to assign
	Matrix& operator=(const Matrix& rhs);

	/* --- FUNCTIONS ---*/
	/// Fills the matrix with given data
	///		@param m[in] matrix represented as an 16-float long 1D array
	///		@param isVec[in] specifies wheather the matrix to be filledi s a 4x4 matrix or 4x1 vector
	void initData(const float *m, bool isVec);

	/// Makes the current matrix an identity matrix
	void makeIdentity();

	/// Function for move constructor
	///		@param rhs[in] reference to the matrix with which data should be swapped 
	void swap(Matrix& rhs);

	/// Outputs the matrix to the console
	void print();
};

#endif