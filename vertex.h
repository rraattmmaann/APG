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

	/// Custom Vertex constructor
	///		@param x[in] vertex x coordinate
	///		@param y[in] vertex y coordinate
	///		@param z[in] vertex z coordinate
	///		@param w[in] vertex w coordinate
	Vertex(float x, float y, float z, float w);

	/// VertexCopy constructor
	Vertex(const Vertex& rhs);

	/* --- OPERATORS ---*/
	/// Operator * for vertex mulitplication with a float value
	///		@param rhs[in] a float value to multiply this vertex(vector) by
	///		@return resulting vertex (vector)
	Vertex operator*(const float& rhs);

	/// Vertex assignment Operator
	///		@param rhs[in] which value to assign
	Vertex& operator=(const Vertex& rhs);

	/* --- FUNCTIONS --- */
	/// Function for move constructor
	///		@param rhs[in] reference to the vertex with which data should be swapped 
	void swap(Vertex& rhs);

	/// Outputs the matrix to the console
	void print();
};
