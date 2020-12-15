#pragma once

#include "vertex.h"

/// Computes dot product of two vectors
///		@param a[in] 1st vector
///		@param b[in] 2nd vector
///		@return scalar value, a dot product of two vectors
inline float dot(const Vertex &a, const Vertex &b);

/// Computes crossproduct of two vertices
///		@param a[in] 1st vector
///		@param b[in] 2nd vector
///		@param res[in] where to store the result
void cross(const Vertex &a, const  Vertex &b, Vertex &res);

/// Normalizes given vector
///		@param a[in] vector to normalize
void normalize(Vertex &a);

/// Returns the distance between 2 vertices (i.e. the length of vector)
///		@param a[in] 1st vector
///		@param b[in] 2nd vector
///		@return scalar value, a distance of given 2 vertices
float vectorLength(const Vertex& a, const Vertex& b);

/// Represents a material of an prpimitive
struct Material {
	float r;
	float g;
	float b;
	float kd;
	float ks;
	float shine;
	float T;
	float ior;

	Material() {
		r = 0.0f;
		g = 0.0f;
		b = 0.0f;
		kd = 0.0f;
		ks = 0.0f;
		shine = 0.0f;
		T = 0.0f;
		ior = 0.0f;
	}

	Material(float _r, float _g, float _b, float _kd, float _ks, float _shine, float _T, float _ior) {
		r = _r;
		g = _g;
		b = _b;
		kd = _kd;
		ks = _ks;
		shine = _shine;
		T = _T;
		ior = _ior;
	}
};

/// Represents a light source in the scene
struct Light {
	Vertex position;
	float r;
	float g;
	float b;

	Light() {}

	Light(Vertex pos, float _r, float _g, float _b) {
		position = pos;
		r = _r;
		g = _g;
		b = _b;
	}

	Light(const Light& rhs) {

		position = rhs.position;
		r = rhs.r;
		g = rhs.g;
		b = rhs.b;
	}

	Light& operator=(const  Light& rhs) {
		Light temp(rhs);
		swap(temp);
		return *this;
	}

	/// Function for move constructor
	void swap(Light& rhs) {
		std::swap(position, rhs.position);
		std::swap(r, rhs.r);
		std::swap(g, rhs.g);
		std::swap(b, rhs.b);
	}
};

/// Represents an emmisive material of a primitive
struct EmissiveMaterial {
	float r;
	float g;
	float b;
	float c0;
	float c1;
	float c2;

	EmissiveMaterial() {}

	EmissiveMaterial(float _r, float _g, float _b, float _c0, float _c1, float _c2) {
		r = _r;
		g = _g;
		b = _b;
		c0 = _c0;
		c1 = _c1;
		c2 = _c2;
	}

};

/// Clss representing an intersection point of a ray and a primitive
struct Intersection {

	/// Point of intersection in global coordinates within the scene
	Vertex position;

	/// Normal of the primitive in the world coordinates in the point of intersection
	Vertex normal;

	/// Distance of the intersection to the camera
	float distance;

	Intersection() : distance(INFINITY) {}

	Intersection(const Intersection& rhs) {
		distance = rhs.distance;
		position = rhs.position;
		normal = rhs.normal;
	}

	Intersection& operator=(const  Intersection& rhs) {
		Intersection temp(rhs);
		swap(temp);
		return *this;
	}

	/// Function for move constructor
	///		@param rhs[in] reference to the intersection with which data should be swapped 
	void swap(Intersection& rhs) {
		std::swap(position, rhs.position);
		std::swap(normal, rhs.normal);
		std::swap(distance, rhs.distance);
	}
};

/// Class represents an enviroment map
struct EnviromentMap {
	int width;
	int height;
	float *texels;

	EnviromentMap() {}

	EnviromentMap(int w, int h, float *tex) {
		width = w;
		height = h;
		texels = tex;
	}
};


/// Class representing a ray casted from the camera into the scene
struct Ray {

	/// The origin of the ray is the camera position
	Vertex origin;

	/// Normalized direction of the ray
	Vertex dir;

	/// Flag indicating refracted ray
	bool refracted = false;

	Ray() {}

	Ray(Vertex _origin, Vertex _dir) {
		origin = _origin;
		dir = _dir;
	}

	Ray(const Ray& rhs) {
		origin = rhs.origin;
		dir = rhs.dir;
		refracted = rhs.refracted;
	}

	Ray& operator=(const  Ray& rhs) {
		Ray temp(rhs);
		swap(temp);
		return *this;
	}

	/// Function for move constructor
	///		@param rhs[in] reference to the ray with which data should be swapped 
	void swap(Ray& rhs) {
		std::swap(origin, rhs.origin);
		std::swap(dir, rhs.dir);
		std::swap(refracted, rhs.refracted);
	}
};

/// Class representing a sphere within the scene
struct Sphere {

	/// Radius of the sphere
	float radius;

	/// Index of the sphere material within the vector of materials
	unsigned int matIdx;

	/// Center of the sphere
	Vertex center;

	Sphere() {}

	Sphere(const Vertex c, float rad, int idx) {
		center = c;
		radius = rad;
		matIdx = idx;
	}

	Sphere(const Sphere& rhs) {
		radius = rhs.radius;
		matIdx = rhs.matIdx;
		center = rhs.center;
	}

	Sphere& operator=(const  Sphere& rhs) {
		Sphere temp(rhs);
		swap(temp);
		return *this;
	}

	/// Function for move constructor
	///		@param rhs[in] reference to the intersection with which data should be swapped 
	void swap(Sphere& rhs) {
		std::swap(radius, rhs.radius);
		std::swap(matIdx, rhs.matIdx);
		std::swap(center, rhs.center);
	}

	/// Computes intersection of this object with given ray
	///		@param ray[in] input ray to check intersection with
	///		@return object of intersection containing information about the intersection, if exists
	Intersection intersects(const Ray &r) {

		float x = center.m_data[0];
		float y = center.m_data[1];
		float z = center.m_data[2];
		Intersection Int;

		Vertex dist(
			r.origin.m_data[0] - x,
			r.origin.m_data[1] - y,
			r.origin.m_data[2] - z,
			1.0f
		);

		float a = dot(r.dir, r.dir);
		float b = 2.0f * dot(dist, r.dir);
		float c = dot(dist, dist) - radius * radius;
		float disc = b * b - 4.0f * a * c;

		float t = -1.0;
		if (disc >= 0) {

			if (r.refracted)
				t = (-b + sqrt(disc)) / (2.0f * a);	// For refracted rays we need the further intersection	
			else 
				t = (-b - sqrt(disc)) / (2.0f * a);	// The closest int. otherwise
				
			if (t < THRESHOLD) return Intersection(); // Eliminates the float inacurracy while computing intersection
			
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];

			Int.normal.m_data[0] = Int.position.m_data[0] - x;
			Int.normal.m_data[1] = Int.position.m_data[1] - y;
			Int.normal.m_data[2] = Int.position.m_data[2] - z;
			normalize(Int.normal);

			Int.distance = t;

			if (dot(r.dir, Int.normal) > 0.0f && !r.refracted) { // Backface culling only for primary and reflected rays
				return Intersection();
			}
		}

		return Int;
	}
};

/// Class representing a polygon within the scene
struct Polygon {
	Vertex a;
	Vertex b;
	Vertex c;
	Vertex normal;
	unsigned int matIdx;
	bool matType;

	Polygon() {}

	Polygon(Vertex _a, Vertex _b, Vertex _c) {
		a = _a;
		b = _b;
		c = _c;
	}

	Polygon(const Polygon& rhs) {
		a = rhs.a;
		b = rhs.b;
		c = rhs.c;
		normal = rhs.normal;
		matIdx = rhs.matIdx;
		matType = rhs.matType;
	}

	Polygon& operator=(const  Polygon& rhs) {
		Polygon temp(rhs);
		swap(temp);
		return *this;
	}

	/// Function for move constructor
	///		@param rhs[in] reference to the intersection with which data should be swapped 
	void swap(Polygon& rhs) {
		std::swap(a, rhs.a);
		std::swap(b, rhs.b);
		std::swap(c, rhs.c);
		std::swap(normal, rhs.normal);
		std::swap(matIdx, rhs.matIdx);
		std::swap(matType, rhs.matType);
	}

	/*
	* https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	*/
	/// Computes intersection of this object with given ray
	///		@param ray[in] input ray to check intersection with
	///		@return object of intersection containing information about the intersection, if exists
	Intersection intersects(Ray &r) {
		Intersection Int;

		if (dot(r.dir, normal) > 0.0f) {
			return Int;
		}

		Vertex edge1, edge2, h, s, q;
		float _a, f, u, v;

		edge1 = b - a;
		edge2 = c - a;
		cross(r.dir, edge2, h);
		_a = dot(edge1, h);
		if (_a > -THRESHOLD && _a < THRESHOLD) {
			return Int;
		}
		f = 1.0f / _a;
		s = r.origin - a;
		u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) {
			return Int;
		}
		cross(s, edge1, q);
		v = f * dot(r.dir, q);
		if (v < 0.0f || u + v > 1.0f) {
			return Int;
		}

		float t = f * dot(edge2, q);
		if (t > THRESHOLD) {
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];
			Int.distance = t;
			Int.normal = normal;			
		}
		return Int;
	}
};

bool operator==(const Polygon& lhs, const Polygon& rhs) {
	if (lhs.a != rhs.a) return false;
	else if (lhs.b != rhs.b) return false;
	else if (lhs.c != rhs.c) return false;
	else if (lhs.normal != rhs.normal) return false;
	return true;
}

inline float dot(const Vertex &a, const Vertex &b)
{
	return a.m_data[0] * b.m_data[0] + a.m_data[1] * b.m_data[1] + a.m_data[2] * b.m_data[2];
}

inline void cross(const Vertex &a, const Vertex &b, Vertex &res) {
	res.m_data[0] = a.m_data[1] * b.m_data[2] - a.m_data[2] * b.m_data[1];
	res.m_data[1] = a.m_data[2] * b.m_data[0] - a.m_data[0] * b.m_data[2];
	res.m_data[2] = a.m_data[0] * b.m_data[1] - a.m_data[1] * b.m_data[0];
}

inline void normalize(Vertex &a) {
	float lenght = sqrt((a.m_data[0] * a.m_data[0]) + (a.m_data[1] * a.m_data[1]) + (a.m_data[2] * a.m_data[2]));
	a.m_data[0] = a.m_data[0] / lenght;
	a.m_data[1] = a.m_data[1] / lenght;
	a.m_data[2] = a.m_data[2] / lenght;
}

float vectorLength(const Vertex& a, const Vertex& b) {
	return sqrt((a.m_data[0] - b.m_data[0]) * (a.m_data[0] - b.m_data[0]) +
		(a.m_data[1] - b.m_data[1]) * (a.m_data[1] - b.m_data[1]) +
		(a.m_data[2] - b.m_data[2]) * (a.m_data[2] - b.m_data[2])
	);
}