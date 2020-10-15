#ifndef MATRIX_H
#define MATRIX_H


class Matrix {
public:
	Matrix(int _height, int _width);
	Matrix() : width(0), height(0), m_data(NULL) {};
	Matrix(const Matrix& rhs);	
	~Matrix();

	//friend std::ostream& operator<<(std::ostream& out, const matrix& m);
	Matrix operator*(const Matrix& rhs);
	Matrix& operator=(const Matrix& rhs);	

	Matrix transpose() const;
	Matrix& multiplySequential(const Matrix& a, const Matrix& b, Matrix& result);
	void makeIdentity();
	void initData(const float *m);

	void swap(Matrix& rhs);
	void allocate();

	unsigned int width;
	unsigned int height;
	float** m_data;
};

/*
std::ostream& operator<<(std::ostream& out, const matrix& m) {

    for (unsigned int i = 0; i < m.height; ++i) {
    	out << "[";
    	for (unsigned int j = 0; j < m.width; ++j) {
    		out << " " << m.m_data[i][j];
    	}
    	out << " ]\n";
    }
    return out;
}*/

#endif