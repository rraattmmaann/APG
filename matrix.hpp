#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
public:
	Matrix(int _height, int _width);
	Matrix() : width(4), height(4) {}; // We are using 4x4 matrices
	Matrix(const Matrix& rhs);	
	~Matrix();

	Matrix operator*(const Matrix& rhs);
	Matrix operator*(const float& rhs);
	Matrix& operator=(const Matrix& rhs);	

	Matrix& multiplySequential(const Matrix& a, const Matrix& b, Matrix& result);
	void makeIdentity();
	void initData(const float *m, bool isVec);
	void print();

	void swap(Matrix& rhs);

	unsigned int width;
	unsigned int height;
	float m_data[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
};

#endif