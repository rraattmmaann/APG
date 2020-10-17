#ifndef MATRIX_H
#define MATRIX_H


class Matrix {
public:
	Matrix(int _height, int _width);
	Matrix() : width(4), height(4) {}; // We are using 4x4 matrices
	Matrix(const Matrix& rhs);	
	~Matrix();

	//friend std::ostream& operator<<(std::ostream& out, const Matrix& m);
	Matrix operator*(const Matrix& rhs);
	Matrix operator*(const float& rhs);
	Matrix& operator=(const Matrix& rhs);	

	Matrix transpose() const;
	Matrix& multiplySequential(const Matrix& a, const Matrix& b, Matrix& result);
	void makeIdentity();
	void initData(const float *m);
	void print();

	void swap(Matrix& rhs);
	void allocate();

	unsigned int width;
	unsigned int height;
	float m_data[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
};


/*std::ostream& operator<<(std::ostream& out, const Matrix& m) {

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