#pragma once
#define _USE_MATH_DEFINES
#include <cmath> 
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "sgl.h"
#include "matrix.hpp"
#include <algorithm>

/// Structure representing a pixel in the canvas
struct colorPixel {
	float r;
	float g;
	float b;
};

struct Material {
	float r;
	float g;
	float b;
	float kd;
	float ks;
	float shine;
	float T;
	float ior;
};

struct Light {
	Vertex position;
	float r;
	float g;
	float b;

	Light(){}

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
	///		@param rhs[in] reference to the intersection with which data should be swapped 
	void swap(Light& rhs) {
		std::swap(position, rhs.position);
		std::swap(r, rhs.r);
		std::swap(g, rhs.g);
		std::swap(b, rhs.b);
	}
};

struct EmissiveMaterial {
	float r;
	float g;
	float b;
	float c0;
	float c1;
	float c2;
};

struct Intersection {
	Vertex position;
	Vertex normal;
	float distance;

	Intersection() : distance(INFINITY){}

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

struct enviromentMap {
	int width;
	int height;
	float *texels;
};

struct Ray {
	Vertex origin;
	Vertex dir;

	Ray() {}

	Ray(const Ray& rhs) {
		origin = rhs.origin;
		dir = rhs.dir;
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
	}
};

float dot(Vertex &a, Vertex &b)
{
	float product = a.m_data[0] * b.m_data[0] + a.m_data[1] * b.m_data[1] + a.m_data[2] * b.m_data[2];
	return product;
}

Vertex cross(Vertex &a, Vertex &b){
	Vertex cross;
	cross.m_data[0] = a.m_data[1] * b.m_data[2] - a.m_data[2] * b.m_data[1];
	cross.m_data[1] = a.m_data[2] * b.m_data[0] - a.m_data[0] * b.m_data[2];
	cross.m_data[2] = a.m_data[0] * b.m_data[1] - a.m_data[1] * b.m_data[0];

	return cross;

}

Vertex normalize(Vertex &a) {
	float lenght = sqrt((a.m_data[0] * a.m_data[0]) + (a.m_data[1] * a.m_data[1]) + (a.m_data[2] * a.m_data[2]));
	a.m_data[0] = a.m_data[0] / lenght;
	a.m_data[1] = a.m_data[1] / lenght;
	a.m_data[2] = a.m_data[2] / lenght;
	return a;
}

Vertex minus(Vertex &a, Vertex &b) {
	Vertex ret;
	ret.m_data[0] = a.m_data[0] - b.m_data[0];
	ret.m_data[1] = a.m_data[1] - b.m_data[1];
	ret.m_data[2] = a.m_data[2] - b.m_data[2];
	return ret;
}

struct Sphere {
	float radius;
	unsigned int matIdx;
	Vertex center;

	Sphere() {}

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

		Vertex normal;

		float a = dot(r.dir, r.dir);
		float b = 2.0 * dot(dist, r.dir);
		float c = dot(dist, dist) - radius * radius;
		float disc = b * b - 4 * a * c;
		
		float t = -1.0;
		if (disc >= 0) {
			t = (-b - sqrt(disc)) / (2.0 * a);
			//mùže být t záporné?
			Int.distance = t;
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];
			
			normal.m_data[0] = Int.position.m_data[0] - x;
			normal.m_data[1] = Int.position.m_data[1] - y;
			normal.m_data[2] = Int.position.m_data[2] - z;
			Int.normal = normalize(normal);
		}

		return Int;
	}
};

struct Polygon {
	Vertex a;
	Vertex b;
	Vertex c;
	Vertex normal;
	unsigned int matIdx;
	bool matType;

	Polygon() {}

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
	Intersection intersects(Ray &r) {
		Intersection Int;
		Vertex vertex0 = a;
		Vertex vertex1 = b;
		Vertex vertex2 = c;
		Vertex edge1, edge2, h, s, q;
		float a, f, u, v;

		edge1 = minus(vertex1, vertex0);
		edge2 = minus(vertex2, vertex0);
		h = cross(r.dir, edge2);
		a = dot(edge1, h);
		if (a > -0.000001 && a < 0.000001) {
			return Int;
		}
		f = 1.0 / a;
		s = minus(r.origin, vertex0);
		u = f * dot(s, h);
		if (u < 0.0 || u > 1.0) {
			return Int;
		}
		q = cross(s, edge1);
		v = f * dot(r.dir, q);
		if (v < 0.0 || u + v > 1.0) {
			return Int;
		}
		
		float t = f * dot(edge2, q);
		if (t > 0.000001) {
			Int.position.m_data[0] = r.origin.m_data[0] + t * r.dir.m_data[0];
			Int.position.m_data[1] = r.origin.m_data[1] + t * r.dir.m_data[1];
			Int.position.m_data[2] = r.origin.m_data[2] + t * r.dir.m_data[2];
			Int.distance = t;
			
			return Int;
		}
		else {


			
		}
		return Int;
	}
};


/// Holds all variables need fot the current context
class Context {
public:
	/* --- VARIABLES --- */
	/// Width and height of the canvas
	int width, height;

	/// Index of the current context
	int index;

	/// Color used for buffer clearing
	colorPixel clearColor;

	/// Color buffer, each pixel is represented as 3 floats, i.e. [ r g b r g b ...]
	float* colorBuffer;
	/// Depth buffer
	float* depthBuffer;

	/// Current viewport matrix
	Matrix viewport;

	/// Current scale determined by the viewport matrix
	float viewportScale;
	
	/// Specifies wheather the depth test is enabled or not
	bool depthTest;

	/// Currently selected matrix mode within this context,
	/// used for switch the access to matrices stack
	enum sglEMatrixMode matrixMode;

	/// Currently rendered primitive type
	enum sglEElementType elementType;

	/// Currently used filling mode of the rendering
	enum sglEAreaMode areaMode;

	/// Currently used color for drawing into the color buffer
	colorPixel drawingColor;

	/// Current size of the rendered point
	float pointSize;

	/// Vector of modelview matrices acting as a stack
	std::vector<Matrix> modelViewMatricesStack;

	/// Vector of projection matrices acting as a stack
	std::vector<Matrix> projectionMatricesStack;

	/// Vector of vertices, fills withing the glBegin() and glEnd() sequence
	std::vector<Vertex> vertexBuffer;

	std::vector<Polygon> polygons;

	std::vector<Light> lights;

	std::vector<EmissiveMaterial> emmisiveMaterials;

	std::vector<enviromentMap> enviromentMaps;

	std::vector<Sphere> spheres;

	std::vector<Material> materials;

	bool addingEmissiveMaterial = false;

	/* --- CONSTRUCTORS & DESTRUCTORS --- */
	/// Default Context constructor
	Context() {}
	
	/// Custom Context constructor
	///		@param _width[in] width of the canvas
	///		@param _height[in] height of the canvas
	///		@param _index[in] the index of current context
	Context(int _width, int _height, int _index) {
		width = _width;
		height = _height;
		index = _index;
		colorBuffer = new float[width * height * 3];
		depthBuffer = new float[width * height];
		viewport.width = 2;
		viewport.height = 4;
		viewport.makeIdentity();		
		depthTest = false;
		clearColor.r = 0;
		clearColor.g = 0;
		clearColor.b = 0;
		drawingColor.r = 1;
		drawingColor.g = 1;
		drawingColor.b = 1;
		clearBuffer((SGL_COLOR_BUFFER_BIT | SGL_DEPTH_BUFFER_BIT));
	}

	/// Default Context destructor
	~Context() {
		try {
			if (colorBuffer) delete[] colorBuffer;
			if (depthBuffer) delete[] depthBuffer;
		}
		catch (const std::exception& e) {
			std::cout << e.what();
		}		

		vertexBuffer.shrink_to_fit();
		projectionMatricesStack.clear();
		modelViewMatricesStack.clear();
		projectionMatricesStack.shrink_to_fit();
		modelViewMatricesStack.shrink_to_fit();
	}
	
	/* --- FUNCTIONS --- */

	/// Computes vertex positions for an Arc
	///		@param x [in] circle center x coordinate
	///		@param y [in] circle center y coordinate
	///		@param z [in] circle center z coordinate
	///		@param a[in] x semiaxis length
	///		@param b[in] y semiaxis length
	void approximationArc(float x, float y, float z, float radius, float from, float to) {

		float x1;
		float y1;
		float x2;
		float y2;
		int steps = 40 * (to - from) / (2 * M_PI);

		if (areaMode == SGL_FILL) {
			sglBegin(SGL_POLYGON);
			sglVertex3f(x, y, z);
		}			
		else
			sglBegin(SGL_LINE_STRIP);

		float alpha = (to - from) / steps;
		float CA = cos(alpha);
		float SA = sin(alpha);

		x1 = radius * cos(from);
		y1 = radius * sin(from);

		sglVertex3f(x + x1, y + y1, z);
		for (int i = 1; i <= steps; i++) {
			x2 = CA * x1 - SA * y1;
			y2 = SA * x1 + CA * y1;
			sglVertex3f(x + x2, y + y2, z);
			x1 = x2;
			y1 = y2;
		}
		sglEnd();
	}

	/// Computes vertex positions for an Ellipse
	///		@param x[in] ellipse center x coordinate
	///		@param y[in] ellipse center y coordinate
	///		@param z[in] ellipse center z coordinate(depth)
	///		@param a[in] x semiaxis length
	///		@param b[in] y semiaxis length
	void approximationEllipse(float x, float y, float z, float a, float b) {

		if (areaMode == SGL_POINT) {
			sglBegin(SGL_POINTS);
			sglVertex3f(x, y, z);
			sglEnd();
		}

		float x2, y2;
		float alpha = 2 * M_PI / 40;
		float x1 = 1;
		float y1 = 0;
		float CA = cos(alpha);
		float SA = sin(alpha);

		if (areaMode == SGL_FILL)
			sglBegin(SGL_POLYGON);
		else
			sglBegin(SGL_LINE_STRIP);

		for (int i = 0; i <= 40; i++) {
			x2 = CA * x1 - SA * y1;
			y2 = SA * x1 + CA * y1;
			sglVertex3f(x + (x2 * a), y + (y2 * b), z);
			x1 = x2;
			y1 = y2;

		}
		sglEnd();
	}

	/// Renders a circle with given radius on defined location
	///		@param xs[in] circle center x coordinate
	///		@param ys[in] circle center y coordinate
	///		@param zs[in] circle center z coordinate
	///		@param r[in] cirlce radius
	void bresenhamCircle(float xs, float ys, float zs, float r) {

		if (areaMode == SGL_POINT) {
			sglBegin(SGL_POINTS);
			sglVertex3f(xs, ys, zs);
			sglEnd();
			return;
		}
		Matrix MV = modelViewMatricesStack.back();
		Matrix P = projectionMatricesStack.back();

		// Compute the transformed circle center
		Vertex stred(xs, ys, zs, 1);
		stred = (P * (MV * stred));
		stred = stred * (1 / stred.m_data[3]);
		int stx = viewport.m_data[0][0] * stred.m_data[0] + viewport.m_data[2][0];
		int sty = viewport.m_data[1][0] * stred.m_data[1] + viewport.m_data[3][0];

		// Compute the transformed radius
		float MVscale = MV.m_data[0][0] * MV.m_data[1][1] - MV.m_data[1][0] * MV.m_data[0][1];
		float Pscale = P.m_data[0][0] * P.m_data[1][1] - P.m_data[1][0] * P.m_data[0][1];
		r *= sqrt(MVscale * Pscale * viewportScale);
		
		if (areaMode == SGL_FILL) {

			int x1;
			int x2;

			int to = (sty + round(r));
			float sqrtContent = (r*r);
			float outerContent = stx + 0.5f;

			for (int count = (sty - r); count <= to; ++count) {
				/*
				// dosazeni do vzorce kruhu a vypocitani xove souradnice bodu - stara verze
				int K = r * r - count * count + 2 * sty * count - sty * sty - stx * stx;
				x2 = (2 * stx + sqrt(4 * stx*stx + 4 * K)) / 2;
				x1 = (2 * stx - sqrt(4 * stx*stx + 4 * K)) / 2;
				*/

				x1 = int(outerContent - sqrt(sqrtContent - ((count - sty)*(count - sty))));
				x2 = int(outerContent + sqrt(sqrtContent - ((count - sty)*(count - sty))));

				bresenhamLine(x1+1, x2, count, count);
			}
		} else {		

			int x, y, p;
			x = 0;
			y = r;
			p = 3 - 2 * r;
			while (x < y) {
				setSymetricalPixels(x, y, stx, sty);
				if (p < 0) {
					p = p + 4 * x + 6;
				}
				else {
					p = p + 4 * (x - y) + 10;
					y = y - 1;
				}
				x = x + 1;
			}
			if (x == y) {
				setSymetricalPixels(x, y, stx, sty);
			}
		}		
	}

	/// Connects 2 given vertices with a line using the Bersenham line algorithm
	///		@param x1[in] vertex No. 1 x coordinate
	///		@param x2[in] vertex No. 2 x coordinate
	///		@param y1[in] vertex No. 1 y coordinate
	///		@param y2[in] vertex No. 2 y coordinate
	void bresenhamLine(int x1, int x2, int y1, int y2) {
		int c0, c1, p;

		if (x2 - x1 <= 0 && y2 - y1 <= 0) {
			int tempX = x1;
			int tempY = y1;

			x1 = x2;
			y1 = y2;

			x2 = tempX;
			y2 = tempY;
		}

		if (abs(y2 - y1) > abs(x2 - x1)) {//svislá
			if (y2 - y1 < 0 && x2 - x1 > 0) {
				int tempX = x1;
				int tempY = y1;

				x1 = x2;
				y1 = y2;

				x2 = tempX;
				y2 = tempY;
			}

			int xDirection = 1;

			if (x2 - x1 <= 0) {
				xDirection = -xDirection;
			}

			c0 = 2 * abs(x2 - x1);
			c1 = c0 - 2 * abs(y2 - y1);
			p = c0 - abs(y2 - y1);

			setPixel(x1, y1);

			for (int i = y1 + 1; i <= y2; i++) {
				if (p < 0) {
					p += c0;
				}
				else {
					p += c1;
					x1 += xDirection;
				}
				setPixel(x1, i);
			}
		}
		else {
			// vodorovná

			if (y2 - y1 >= 0 && x2 - x1 <= 0) {
				int tempX = x1;
				int tempY = y1;

				x1 = x2;
				y1 = y2;

				x2 = tempX;
				y2 = tempY;
			}

			setPixel(x1, y1);
			int yDirection = 1;

			if (y2 - y1 <= 0) {
				yDirection = -yDirection;
			}

			c0 = 2 * abs(y1 - y2);
			c1 = c0 - 2 * abs(x2 - x1);
			p = c0 - abs(x2 - x1);

			for (int i = x1 + 1; i < x2; i++) {
				if (p < 0) {
					p += c0;
				}
				else {
					p += c1;
					y1 += yDirection;
				}
				setPixel(i, y1);
			}
		}
	}

	/// Takes all vertices from the vertex buffer and handles their rendering
	void drawPoints() {

		Matrix PVM = projectionMatricesStack.back() * modelViewMatricesStack.back();

		for (Vertex vert : vertexBuffer) {

			Vertex res = PVM * vert;
			res = res * (1 / res.m_data[3]);

			int x = static_cast<int>(viewport.m_data[0][0] * res.m_data[0] + viewport.m_data[2][0]);
			int y = static_cast<int>(viewport.m_data[1][0] * res.m_data[1] + viewport.m_data[3][0]);

			int currPointSize = pointSize;

			if (currPointSize % 2 == 0) {
				x = x - currPointSize * 0.5 - 1;
				y = y - currPointSize * 0.5 - 1;
			}
			else {
				x = x - (currPointSize - 1) * 0.5;
				y = y - (currPointSize - 1) * 0.5;
			}

			for (int a = 0; a < currPointSize; a++) {
				for (int b = 0; b < currPointSize; b++) {
					setPixel(x + a, y + b);
				}
			}
		}
	}

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	void drawLines() {

		if (vertexBuffer.size() % 2 > 0) return;

		Matrix PVM = projectionMatricesStack.back() * modelViewMatricesStack.back();

		for (unsigned int i = 0; i < vertexBuffer.size(); i += 2) {
			Vertex v1 = vertexBuffer[i];
			Vertex v2 = vertexBuffer[i + 1];
			v1 = PVM * v1;
			v1 = v1 * (1 / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1 / v2.m_data[3]);

			bresenhamLine(
				viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0],
				viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]
			);
		}
	}

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	void drawLineStrip() {

		Matrix PVM = projectionMatricesStack.back() * modelViewMatricesStack.back();

		for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {

			Vertex v1 = vertexBuffer[i];
			Vertex v2 = vertexBuffer[i + 1];
			v1 = PVM * v1;
			v1 = v1 * (1 / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1 / v2.m_data[3]);

			bresenhamLine(
				viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0],
				viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]
			);
		}
	}

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	/// and connects the first and the last vertex
	void drawLineLoop() {

		if (vertexBuffer.size() == 0) return;

		Matrix PVM = projectionMatricesStack.back() * modelViewMatricesStack.back();

		Vertex v = vertexBuffer[0]; // Because of this access to the vertex buffer we need to check if it contains at least 1 alement
		v = PVM * v;
		v = v * (1 / v.m_data[3]);

		int startx = viewport.m_data[0][0] * v.m_data[0] + viewport.m_data[2][0];
		int starty = viewport.m_data[1][0] * v.m_data[1] + viewport.m_data[3][0];

		int length = 0;

		for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {
			Vertex v1 = vertexBuffer[i];
			Vertex v2 = vertexBuffer[i + 1];
			v1 = PVM * v1;
			v1 = v1 * (1 / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1 / v2.m_data[3]);

			bresenhamLine(
				viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0],
				viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0],
				viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]
			);
			length++;
		}

		Vertex v2 = vertexBuffer[length];
		v2 = PVM * v2;
		v2 = v2 * (1 / v2.m_data[3]);

		int endx = viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0];
		int endy = viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0];

		bresenhamLine(endx, startx, endy, starty);
	}

	/// Finds intersection in the x axis in height of y on line between vertices a and b
	///		@param y[in] current height of the scanline
	///		@param a[in] bottom end of the edge
	///		@param b[in] upper end of the edge
	///		@param intersections[in] vector of x-coordinates of all found intersections
	inline void getIntersections(int y, Vertex &a, Vertex &b, std::vector<int> &intersections) {
		int x1 = a.m_data[0];
		int y1 = a.m_data[1];
		int x2 = b.m_data[0];		
		int y2 = b.m_data[1];

		if (x2-x1 == 0) {
			intersections.emplace_back(x1);
			return;
		}

		float slope = (float)(y2 - y1) / (float)(x2 - x1);
		
		float shift = y1 - slope * x1;

		// y = slope * x + shift
		// dosadime za y promenou y xd
		int x = (y - shift) / slope;
		intersections.emplace_back(x);
	}

	/// Finds intersection in the x axis in height of y on line between vertices a and b
	/// Computes the z value of the intersection as well
	///		@param y[in] current height of the scanline
	///		@param a[in] bottom end of the edge
	///		@param b[in] upper end of the edge
	///		@param intersections[in] vector of couple (x-coordinate, depth) of all found intersections
	inline void getIntersectionsDepth(int y, Vertex &a, Vertex &b, std::vector<std::tuple<int, float>> &intersections) {
		int x1 = a.m_data[0];
		int y1 = a.m_data[1];
		int x2 = b.m_data[0];
		int y2 = b.m_data[1];
		float z1 = a.m_data[2];
		float z2 = b.m_data[2];
		float z = (z1*(y2 - y) + z2 * (y - y1)) / (y2 - y1);

		if (x2 - x1 == 0) {
			intersections.emplace_back(std::make_tuple(x1, z));
			return;
		}

		float slope = (float)(y2 - y1) / (float)(x2 - x1);

		float shift = y1 - slope * x1;

		// y = slope * x + shift
		// dosadime za y promenou y xd
		int x = (y - shift) / slope;
		intersections.emplace_back(std::make_tuple(x, z));
	}

	/// Draws polygon, filling depends on the currently set area mode
	///  - SGL_FILL - polygon is filled using scanline filling algorithm
	///  - SGL_LINE - only the polygon edges are rendered using bresenham in line loop
	void drawPolygon() {
		if (areaMode == SGL_FILL) {

			std::vector<Vertex> screenSpaceVertices;
			int yMax = 0;
			int yMin = height;

			std::vector<std::tuple<int, float>> intersectionsDepth;
			std::vector<int> intersections;

			Matrix PVM = projectionMatricesStack.back() * modelViewMatricesStack.back();
			int length = 0;
			
			// Get all screenspace vertices of the polygon first
			for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {
				Vertex v1 = vertexBuffer[i];
				Vertex v2 = vertexBuffer[i + 1];
				v1 = PVM * v1;
				v1 = v1 * (1 / v1.m_data[3]);
				v2 = PVM * v2;
				v2 = v2 * (1 / v2.m_data[3]);

				int y1 = viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0];
				int y2 = viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0];
				
				screenSpaceVertices.emplace_back(Vertex(
					viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0],
					y1,
					viewport.m_data[0][1] * v1.m_data[2] + viewport.m_data[1][1],
					1)
				);
				screenSpaceVertices.emplace_back(Vertex(
					viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0],
						y2, 
					viewport.m_data[0][1] * v2.m_data[2] + viewport.m_data[1][1],
						1)
				);
				length++;

				yMax = std::max(std::max(y1, y2), yMax);
				yMin = std::min(std::min(y1, y2), yMin);
			}

			// Add first and the last vertex into the list
			Vertex v = vertexBuffer[0];
			v = PVM * v;
			v = v * (1 / v.m_data[3]);
			Vertex v2 = vertexBuffer[length];
			v2 = PVM * v2;
			v2 = v2 * (1 / v2.m_data[3]);

			int startx = viewport.m_data[0][0] * v.m_data[0] + viewport.m_data[2][0];
			int starty = viewport.m_data[1][0] * v.m_data[1] + viewport.m_data[3][0];
			float startz = viewport.m_data[0][1] * v.m_data[2] + viewport.m_data[1][1];		

			int endx = viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0];
			int endy = viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0];
			float endz = viewport.m_data[0][1] * v2.m_data[2] + viewport.m_data[1][1];

			screenSpaceVertices.emplace_back(Vertex(endx, endy, endz, 1));
			screenSpaceVertices.emplace_back(Vertex(startx, starty, startz, 1));

			// Run the scanline algorithm on the screenspace vertices
			for (int y = yMax; y > yMin; --y) {
				
				// Find all intersections of the current scanline with the edges
				// POZNAMKA: toto lze zrychlit predpocitanim hran a tedy i jejich smernic(kazda hrana je dana prinkou s predpisem y = ax + b)
				//           protoze vypocet se stale opakuje pro kazdou uroven scanline
				for (unsigned int i = 0; i < screenSpaceVertices.size(); i+=2) {

					Vertex a = screenSpaceVertices[i];
					Vertex b = screenSpaceVertices[i + 1];

					int y1 = a.m_data[1];
					int y2 = b.m_data[1];

					if (y1 < y2) {
						int tmp = y1;
						y1 = y2;
						y2 = tmp;
					}

					if (y1 == y2 || (y > y1 && y > y2) || (y < y1 && y <= y2)) continue;

					if (y1 > y2)
						(depthTest) ? getIntersectionsDepth(y, a, b, intersectionsDepth) : getIntersections(y, a, b, intersections);
					else
						(depthTest) ? getIntersectionsDepth(y, b, a, intersectionsDepth) : getIntersections(y, b, a, intersections);
				}
						
				// Check for intersections
				if (intersections.size() > 1 || intersectionsDepth.size() > 1) {
					// Render with Z buffer
					if (depthTest) {
						std::sort(intersectionsDepth.begin(), intersectionsDepth.end());
						for (unsigned int i = 0; i < intersectionsDepth.size() - 1; i += 2) {

							int x1 = std::get<0>(intersectionsDepth[i]);
							int x2 = std::get<0>(intersectionsDepth[i + 1]);
							float z1 = std::get<1>(intersectionsDepth[i]);
							float z2 = std::get<1>(intersectionsDepth[i + 1]);

							float temp = (x2*z1 - z2 * x1);

							for (int j = x1 + 1; j <= x2; ++j) {
								float z = 1/((z2*j - z1*j + temp) / (x2 - x1));
								int idx = y * width + j;
								if (depthBuffer[idx] > z) {
									depthBuffer[idx] = z;
									setPixel(j, y);									
								}
							}
						}
						intersectionsDepth.clear();
					}
					else { // Render without Z buffer
						std::sort(intersections.begin(), intersections.end());
						for (unsigned int i = 0; i < intersections.size() - 1; i += 2) {

							int x1 = intersections[i];
							int x2 = intersections[i + 1];

							for (int j = x1 + 1; j <= x2; ++j) {
								setPixel(j, y);
							}
						}
						intersections.clear();
					}					
				}
			}
		}
		else { // areaMode = SGL_LINE
			drawLineLoop();
		}
	}

	/// Draws area light in the scene
	void drawAreaLight() {
		// TODO
	}
	

	/// Places currently used color to the specifed position in the color buffer
	///		@param x[in] pixel x coordinate
	///		@param y[in] pixel y coordinate
	inline void setPixel(int x, int y) {
		if (y >= height || y < 0 || x < 0 || x > width) {
			return;
		}

		int position = x + y * width;
		position *= 3;
		colorBuffer[position] = drawingColor.r;
		colorBuffer[position + 1] = drawingColor.g;
		colorBuffer[position + 2] = drawingColor.b;
	}

	/// Places given color to the specifed position in the color buffer
	///		@param x[in] pixel x coordinate
	///		@param y[in] pixel y coordinate
	inline void setPixel(int x, int y, float r, float g, float b) {
		if (y >= height || y < 0 || x < 0 || x > width) {
			return;
		}

		int position = x + y * width;
		position *= 3;
		colorBuffer[position] = r;
		colorBuffer[position + 1] = g;
		colorBuffer[position + 2] = b;
	}

	/// Places currently used color to 8 specifed positions in the color buffer
	/// while rendering a circle
	///		@param x[in] current x axis offset
	///		@param y[in] current y axis offset
	///		@param xs[in] circle center x coordinate
	///		@param ys[in] circle center y coordinate
	void setSymetricalPixels(int x, int y, int xs, int ys) {
		setPixel(xs + x, ys + y);
		setPixel(xs + y, ys + x);
		setPixel(xs + y, ys - x);
		setPixel(xs + x, ys - y);
		setPixel(xs - x, ys - y);
		setPixel(xs - y, ys - x);
		setPixel(xs - y, ys + x);
		setPixel(xs - x, ys + y);
	}
	
	/// Sets current viewport matrix
	///		@param x[in] viewport origin x coordinate(with respect to the canvas)
	///		@param y[in] viewport origin y coordinate(with respect to the canvas)
	///		@param width[in] viewport width in pixels
	///		@param height[in] viewport height in pixels
	void setViewport(float x, float y, float width, float height) {
		viewport.m_data[0][0] = width / 2.0f;
		viewport.m_data[1][0] = height / 2.0f;
		viewport.m_data[2][0] = x + width / 2.0f;
		viewport.m_data[3][0] = y + height / 2.0f;

		viewportScale = (width*height) / 4;
	}

	/// Clear the buffer specified by the input value, returns error otherwise
	///		@param what[in] type of buffer to clear
	///		@return SGL_NO_ERROR on successful clear, SGL_INVALID_VALUE otherwise
	sglEErrorCode clearBuffer(unsigned what) {
		if (what == (SGL_COLOR_BUFFER_BIT | SGL_DEPTH_BUFFER_BIT)) {
			for (int i = 0; i < width * height; i += 3) {
				colorBuffer[i] = clearColor.r;
				colorBuffer[i + 1] = clearColor.g;
				colorBuffer[i + 2] = clearColor.b;
				depthBuffer[i] = depthBuffer[i+1] = depthBuffer[i+2] = INFINITY;
			}
		} else if (what == SGL_COLOR_BUFFER_BIT) {
			for (int i = 0; i < width * height; i += 3) {
				colorBuffer[i] = clearColor.r;
				colorBuffer[i + 1] = clearColor.g;
				colorBuffer[i + 2] = clearColor.b;
			}
		}
		else if (what == SGL_DEPTH_BUFFER_BIT) {
			for (int i = 0; i < width * height; ++i) {
				depthBuffer[i] = INFINITY;
			}
		}
		else {
			return SGL_INVALID_VALUE;
		}
		return SGL_NO_ERROR;
	}

	float Q_rsqrt(float number) {
		long i;
		float x2, y;
		const float threehalfs = 1.5F;

		x2 = number * 0.5F;
		y = number;
		i = *(long *)&y;                       // evil floating point bit level hacking
		i = 0x5f3759df - (i >> 1);               // what the fuck? 
		y = *(float *)&i;
		y = y * (threehalfs - (x2 * y * y));   // 1st iteration
	//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

		return y;
	}

	void Normalize(Vertex &v) {
		float length = Q_rsqrt(
			(v.m_data[0] * v.m_data[0]) +
			(v.m_data[1] * v.m_data[1]) +
			(v.m_data[2] * v.m_data[2]));
		v.m_data[0] /= length;
		v.m_data[1] /= length;
		v.m_data[2] /= length;
	}


	Vertex phong(Ray ray, Intersection intersection, Material mat) {
		Vertex ret;
		float kd, ks, shine, r,g,b;
		kd = mat.kd;
		ks = mat.ks;
		shine = mat.shine;
		r = mat.r;
		g = mat.g;
		b = mat.b;

		for (unsigned int i = 0; i < lights.size(); i++) {
			Light l = lights[i];
			Vertex v;
			Vertex distance = minus(l.position, intersection.position) * minus(l.position, intersection.position);
			float distanceF = sqrt(distance.m_data[0] + distance.m_data[1] + distance.m_data[2]);
			
			Vertex L = normalize(minus( l.position, intersection.position));
			Vertex minusL = L*-1;
			Vertex R = minus(minusL, intersection.normal * dot(minusL, intersection.normal) * 2);
			Vertex V = ray.dir * -1;

			Vertex matColor;
			matColor.m_data[0] = mat.r;
			matColor.m_data[1] = mat.g;
			matColor.m_data[2] = mat.b;

			Vertex lightColor;
			lightColor.m_data[0] = l.r;
			lightColor.m_data[1] = l.g;
			lightColor.m_data[2] = l.b;

			//v += matColor * lightColor; nechteji po nas
			v += matColor  * std::max(0.0f, dot(intersection.normal, L)) * mat.kd;
			v += std::pow(std::max(0.0f, dot(R, V)), mat.shine) * mat.ks;

			//v = v * (1.0/distanceF); nezname intenzitu svetla, takze nedavame
			ret += v;
		}
		
		return ret;
	}

	void startRt(){

		Matrix MV = modelViewMatricesStack.back();
		Matrix P = projectionMatricesStack.back();
		Matrix invMV = MV.inverse();
		Matrix invP = P.inverse();
		Matrix Vp;
		Vp.m_data[0][0] = viewport.m_data[0][0];
		Vp.m_data[0][3] = viewport.m_data[2][0];
		Vp.m_data[1][1] = viewport.m_data[1][0];
		Vp.m_data[1][3] = viewport.m_data[3][0];
		Vp.m_data[2][2] = 1; // 22
		Vp.m_data[2][3] = 0; // 23
		Matrix invVp = Vp.inverse();
		//invVp.m_data[2][2] = 1;
		//invVp.m_data[2][3] = 0;
		Matrix MVP = invMV * P.inverse() * invVp;

		

		// transformovat body do world coords TODO
		// vypocitat normaly polygonu TODO
		/*for (unsigned int i = 0; i < polygons.size(); ++i) {
			polygons[i].a = MV * polygons[i].a;
			polygons[i].b = MV * polygons[i].b;
			polygons[i].c = MV * polygons[i].c;
			Vertex edge1 = minus(polygons[i].b, polygons[i].a);
			Vertex edge2 = minus(polygons[i].c, polygons[i].a);
			polygons[i].normal = cross(edge1, edge2);
			Normalize(polygons[i].normal);
		}

		for (unsigned int i = 0; i < spheres.size(); ++i) {
			spheres[i].center = MV * spheres[i].center;
			float MVscale = MV.m_data[0][0] * MV.m_data[1][1] - MV.m_data[1][0] * MV.m_data[0][1];
			spheres[i].radius *= Q_rsqrt(MVscale);
		}*/

		

		// urcit paprsky TODO
		Vertex rayOrigin = invMV * Vertex(0, 0, 0, 1);
		
		// iff prunik -> vypocit osvetleni ze vsech point lights a nastavit barvu
		// neni prunik -> barva pozadí
		
		
		
		// projit vsechny pixely a vrhat paprsky
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				
				// urcit paprsky TODO
				Ray r;
				r.origin = rayOrigin;
				r.dir = normalize(minus(MVP * Vertex(x, y, -1, 1), rayOrigin));

				Intersection bestInt;
				Polygon bestPolygon;
				Sphere bestSphere;
				bool polygonWins = true;

				// kazdy paprsek zkusit prunik se scenou
				for (unsigned int i = 0; i < polygons.size(); ++i) {
					Intersection Int = polygons[i].intersects(r);
					if (Int.distance < bestInt.distance) {
						bestInt = Int;
						bestPolygon = polygons[i];
					}
				}

				for (unsigned int i = 0; i < spheres.size(); ++i) {
					Intersection Int = spheres[i].intersects(r);
					if (Int.distance < bestInt.distance) {
						bestInt = Int;
						bestSphere = spheres[i];
						polygonWins = false;
					}
				}

				if (bestInt.distance < INFINITY) {
					// nalezli jsme prusecik nejblize kamery s danym primitivem
					// pripocitat svetlo a nakreslit do FB
					if (polygonWins) {
						// polygon
						if (!bestPolygon.matType) {
							Vertex v = phong(r, bestInt, materials[bestPolygon.matIdx]);

							setPixel(x, y,
								v.m_data[0],
								v.m_data[1],
								v.m_data[2]
							);
						}
						else {
							//TODO
							Vertex v = phong(r, bestInt, materials[bestPolygon.matIdx]);

							setPixel(x, y,
								v.m_data[0], 
								v.m_data[1],
								v.m_data[2]
							);
							
						}
					}
					else {
						// sphere

						Vertex v = phong(r, bestInt, materials[bestSphere.matIdx]);

						setPixel(x, y,
							v.m_data[0] ,
							v.m_data[1] ,
							v.m_data[2]
						);
					}					
				}
				else {
					// nastavit do FB pozadi
					setPixel(x, y,
						clearColor.r,
						clearColor.g,
						clearColor.b
					);
				}
			}
		}

	}

	void storePolygons() {
		
		Polygon p;
		p.a = vertexBuffer[0];
		p.b = vertexBuffer[1];
		p.c = vertexBuffer[2];
		Vertex edge1 = minus(p.b, p.a);
		Vertex edge2 = minus(p.c, p.a);
		p.normal = cross(edge1, edge2);
		p.normal = normalize(p.normal);
		p.matIdx = (addingEmissiveMaterial) ? emmisiveMaterials.size() - 1 : materials.size() - 1;
		p.matType = addingEmissiveMaterial;
		polygons.push_back(p);
	}
};