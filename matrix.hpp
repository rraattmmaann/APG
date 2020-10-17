#ifndef MATRIX_H
#define MATRIX_H

/// Class representing a 4x4 flat matrix, can be used as a representation of points
class Matrix {
public:
	Matrix(int _height, int _width);
	Matrix() : width(4), height(4) {}; // We are using 4x4 matrices
	Matrix(const Matrix& rhs);	
	~Matrix();

	Matrix operator*(const Matrix& rhs);
	Matrix operator*(const float& rhs);
	Matrix& operator=(const Matrix& rhs);
	void swap(Matrix& rhs);

	/// Multiplies given matrices and returns the resulting matrix
	///		@param a[in] lhs matrix
	///		@param b[in] rhs matrix
	///		@param result[in] result matrix
	///		@return reference to the result matrix
	Matrix& multiplySequential(const Matrix& a, const Matrix& b, Matrix& result);

	/// Makes the current matrix an identity matrix
	void makeIdentity();

	/// Fills the matrix with given data
	///		@param m[in] matrix represented as an 16-float long 1D array
	///		@param isVec[in] specifies wheather the matrix to be filledi s a 4x4 matrix or 4x1 vector
	void initData(const float *m, bool isVec);

	/// Outputs the matrix to the console
	void print();

	/// Width of the matrix
	unsigned int width;

	/// Height of the matrix
	unsigned int height;

	/// The matrix itself
	float m_data[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
};

#endif