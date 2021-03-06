#pragma once
#include <iostream>

/// Class representing a 4D vertex, or a 4D vector
class Vertex {
public:
	/* --- VARIABLES --- */
	/// Vertex data
	float m_data[4] = { 0,0,0,0 };

	/* --- CONSTRUCTORS & DESTRUCTORS --- */
	/// Default Vertex constructor
	Vertex() {};

	/// Default Vertex constructor
	~Vertex() = default;

	/// Custom Vertex constructor
	///		@param x[in] vertex x coordinate
	///		@param y[in] vertex y coordinate
	///		@param z[in] vertex z coordinate
	///		@param w[in] vertex w coordinate
	Vertex(float x, float y, float z, float w) {
		m_data[0] = x;
		m_data[1] = y;
		m_data[2] = z;
		m_data[3] = w;
	}

	/// VertexCopy constructor
	Vertex(const Vertex& rhs) {
		std::copy(rhs.m_data, rhs.m_data + 4, m_data);
	}

	/* --- OPERATORS ---*/
	/// Operator * for vertex mulitplication with a float value
	///		@param rhs[in] a float value to multiply this vertex(vector) by
	///		@return resulting vertex (vector)
	Vertex operator*(const float& rhs) {
		Vertex result;
		for (unsigned int i = 0; i < 4; i++) {
			result.m_data[i] = m_data[i] * rhs;
		}
		return result;
	}

	Vertex operator*(const Vertex& rhs) {
		Vertex result;
		for (unsigned int i = 0; i < 4; i++) {
			result.m_data[i] = m_data[i] * rhs.m_data[i];
		}
		return result;
	}

	Vertex& operator+=(const Vertex& rhs) {

		m_data[0] += rhs.m_data[0];
		m_data[1] += rhs.m_data[1];
		m_data[2] += rhs.m_data[2];
		m_data[3] += rhs.m_data[3];
		return *this;
	}

	Vertex& operator+=(const float& rhs) {

		m_data[0] += rhs;
		m_data[1] += rhs;
		m_data[2] += rhs;
		m_data[3] += rhs;
		return *this;
	}

	/// Vertex assignment Operator
	///		@param rhs[in] which value to assign
	Vertex& operator=(const Vertex& rhs) {
		Vertex temp(rhs);
		swap(temp);
		return *this;
	}

	/// Adds this vector from vector rhs
	///		@param rhs[in] 2nd vector
	///		@return result of addition of this vector and rhs
	Vertex operator+(const Vertex& rhs) {
		return Vertex(m_data[0] + rhs.m_data[0], m_data[1] + rhs.m_data[1], m_data[2] + rhs.m_data[2], 1);
	}

	/// Substracts this vector from vector rhs
	///		@param rhs[in] 2nd vector
	///		@return result of substraction of this vector and rhs
	Vertex operator-(const Vertex& rhs) {
		return Vertex(m_data[0] - rhs.m_data[0], m_data[1] - rhs.m_data[1], m_data[2] - rhs.m_data[2], 1);
	}

	

	/* --- FUNCTIONS --- */
	/// Function for move constructor
	///		@param rhs[in] reference to the vertex with which data should be swapped 
	void swap(Vertex& rhs) {
		std::swap(m_data, rhs.m_data);
	}

	/// Outputs the matrix to the console
	void print() {
		std::cout << "[";
		for (unsigned int j = 0; j < 4; ++j) {
			std::cout << " " << m_data[j];
		}
		std::cout << " ]\n";
	}
};

bool operator==(const Vertex& lhs, const Vertex& rhs) {
	if (lhs.m_data[0] != rhs.m_data[0]) return false;
	else if (lhs.m_data[1] != rhs.m_data[1]) return false;
	else if (lhs.m_data[2] != rhs.m_data[2]) return false;
	return true;
}

bool operator!=(const Vertex& lhs, const Vertex& rhs) {
	return !(rhs == lhs);
}
