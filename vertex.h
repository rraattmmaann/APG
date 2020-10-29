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

	/// Vertex assignment Operator
	///		@param rhs[in] which value to assign
	Vertex& operator=(const Vertex& rhs) {
		Vertex temp(rhs);
		swap(temp);
		return *this;
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
