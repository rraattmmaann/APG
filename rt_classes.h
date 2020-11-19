#pragma once

#include "vertex.h"

/// Computation of sqrt of a floating point number
/// ! This function was taken fom the Quake 3 game implementation !
/// We don't take any credits for this
///		@param number[in] number to compute sqrt from
float Q_rsqrt(float number);

/// Computes dot product of two vectors
///		@param a[in] 1st vector
///		@param b[in] 2nd vector
///		@return scalar value, a dot product of two vectors
inline float dot(Vertex &a, Vertex &b);

/// Computes crossproduct of two vertices
///		@param a[in] 1st vector
///		@param b[in] 2nd vector
///		@param res[in] where to store the result
void cross(Vertex a, Vertex b, Vertex &res);

/// Normalizes given vector
///		@param a[in] vector to normalize
void normalize(Vertex &a);

/// Represents a material of an primitive
struct Material {
	float r;
	float g;
	float b;
	float kd;
	float ks;
	float shine;
	float T;
	float ior;

	Material() {}

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

	int depth;

	Ray() {
		depth = 8;
	}

	Ray(const Ray& rhs) {
		origin = rhs.origin;
		dir = rhs.dir;
		depth = rhs.depth;
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
		std::swap(depth, rhs.depth);
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

	Sphere(Vertex c, float rad, float idx) {
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
	Intersection intersects(Ray &r) {

		float x = center.m_data[0];
		float y = center.m_data[1];
		float z = center.m_data[2];
		Intersection Int;
		//Vertex dist = r.origin - center;
		Vertex dist;
		dist.m_data[0] = r.origin.m_data[0] - x;
		dist.m_data[1] = r.origin.m_data[1] - y;
		dist.m_data[2] = r.origin.m_data[2] - z;

		float a = dot(r.dir, r.dir);
		float b = 2.0 * dot(dist, r.dir);
		float c = dot(dist, dist) - radius * radius;
		float disc = b * b - 4 * a * c;

		float t = -1.0;
		if (disc >= 0) {
			t = (-b - sqrt(disc)) / (2.0 * a);
			//mùe bıt t záporné?
			Int.distance = t;
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];

			Int.normal.m_data[0] = Int.position.m_data[0] - x;
			Int.normal.m_data[1] = Int.position.m_data[1] - y;
			Int.normal.m_data[2] = Int.position.m_data[2] - z;
			normalize(Int.normal);
			

			if (dot(r.dir, Int.normal) > 0.0) {
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
		//Vertex vertex0 = a;
		//Vertex vertex1 = b;
		//Vertex vertex2 = c;

		if (dot(r.dir, normal) > 0.0) {
			return Int;
		}

		Vertex edge1, edge2, h, s, q;
		float _a, f, u, v;

		edge1 = b - a;
		edge2 = c - a;
		cross(r.dir, edge2, h);
		_a = dot(edge1, h);
		if (_a > -0.01 && _a < 0.01) {
			return Int;
		}
		f = 1.0 / _a;
		s = r.origin - a;
		u = f * dot(s, h);
		if (u < 0.0 || u > 1.0) {
			return Int;
		}
		cross(s, edge1, q);
		v = f * dot(r.dir, q);
		if (v < 0.0 || u + v > 1.0) {
			return Int;
		}

		float t = f * dot(edge2, q);
		if (t > 0.001) {
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];
			Int.distance = t;
			Int.normal = normal;			
		}
		return Int;
	}
};

float Q_rsqrt(float number) {
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y = number;
	i = *(long *)&y;                       // evil floating point bit level hacking
	i = 0x5f3759df - (i >> 1);             // what the fuck? 
	y = *(float *)&i;
	y = y * (threehalfs - (x2 * y * y));   // 1st iteration
//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}

inline float dot(Vertex &a, Vertex &b)
{
	return a.m_data[0] * b.m_data[0] + a.m_data[1] * b.m_data[1] + a.m_data[2] * b.m_data[2];
}

inline void cross(Vertex a, Vertex b, Vertex &res) {
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