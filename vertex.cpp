#include "vertex.h"

Vertex::Vertex(float x, float y, float z, float w) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

Vertex::Vertex(const Vertex& rhs) {
	std::copy(rhs.m_data, rhs.m_data + 4, m_data);
}

Vertex Vertex::operator*(const float& rhs) {
	Vertex result;
	for (unsigned int i = 0; i < 4; i++) {
		result.m_data[i] = m_data[i] * rhs;
	}
	return result;
}

Vertex& Vertex::operator=(const Vertex& rhs) {
	Vertex temp(rhs);
	swap(temp);
	return *this;
}

void Vertex::swap(Vertex& rhs) {
	std::swap(m_data, rhs.m_data);
}

void Vertex::print() {
	std::cout << "[";
	for (unsigned int j = 0; j < 4; ++j) {
		std::cout << " " << m_data[j];
	}
	std::cout << " ]\n";
}