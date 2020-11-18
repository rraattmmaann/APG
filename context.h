#pragma once
#define _USE_MATH_DEFINES
#include <cmath> 
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <thread>

#include "sgl.h"
#include "matrix.hpp"
#include "rt_classes.h"

/// Structure representing a pixel in the canvas
struct colorPixel {
	float r;
	float g;
	float b;
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

	std::vector<EnviromentMap> enviromentMaps;

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
		r *= Q_rsqrt(MVscale * Pscale * viewportScale);
		
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

				x1 = int(outerContent - Q_rsqrt(sqrtContent - ((count - sty)*(count - sty))));
				x2 = int(outerContent + Q_rsqrt(sqrtContent - ((count - sty)*(count - sty))));

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

		if (abs(y2 - y1) > abs(x2 - x1)) {//svisl�
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
			// vodorovn�

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

	/// Phong illumination model
	///		@param ray[in] ray casted from the viewer(camera) to the primitive
	///		@param intersection[in] the point of intersection of the ray with the primitive
	///		@param mat[in] material used for the lightning model
	///		@return vertex cointaining the resulting color of the pixel (r,g,b,1)
	Vertex phong(Ray &ray, Intersection &intersection, Material &mat, float xPos) {
		Vertex ret;
		Vertex matColor(mat.r, mat.g, mat.b, 1);
		matColor = matColor * mat.kd;

		

		for (unsigned int i = 0; i < lights.size(); i++) {
			Light l = lights[i];
			Vertex toLight = minus(intersection.position, l.position);
	
			normalize(toLight);
			Ray r;
			r.dir = toLight*-1;
			r.origin = intersection.position;
			bool ok = true;

			for (unsigned int j = 0; j < polygons.size(); ++j) {
				Intersection Int = polygons[j].intersects(r);
				if (Int.distance != INFINITY) {
					ok = false;
					break;
				}
			}
			if (!ok) {
				continue;
			}
			for (unsigned int j = 0; j < spheres.size(); ++j) {
				if (spheres[j].center.m_data[0] == xPos) {
					continue;
				}
				Intersection Int = spheres[j].intersects(r);
				if (Int.distance != INFINITY) {
					ok = false;
					break;
				}
			}
			if (!ok) {
				continue;
			}
			
			Vertex light(l.r, l.g, l.b, 1);

			Vertex L = minus( l.position, intersection.position);
			normalize(L);
			Vertex minusL = L*-1;
			Vertex R = minus(minusL, intersection.normal * dot(minusL, intersection.normal) * 2);

			Vertex v = matColor * light * std::max(0.0f, dot(intersection.normal, L));
			v += std::pow(std::max(0.0f, dot(R, ray.dir)), mat.shine) * mat.ks;

			ret += v;
		}
		
		return ret;
	}

	/// Starts ray tracing in the current scene
	void startRt(){

		// get all matrices ready
		Matrix invMV = modelViewMatricesStack.back().inverse();
		Matrix Vp;
		Vp.m_data[0][0] = viewport.m_data[0][0];
		Vp.m_data[0][3] = viewport.m_data[2][0];
		Vp.m_data[1][1] = viewport.m_data[1][0];
		Vp.m_data[1][3] = viewport.m_data[3][0];
		Matrix MVP = invMV * projectionMatricesStack.back().inverse() * Vp.inverse();

		// define the origin of all the rays
		Ray r;
		r.origin = invMV * Vertex(0, 0, 0, 1);
		
		// go through all of the pixels and cast rays from the camera through each pixel
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				
				// direction of current ray				
				r.dir = minus(MVP * Vertex(x , y, -1, 1), r.origin);
				normalize(r.dir);

				Intersection bestInt;
				Polygon bestPolygon;
				Sphere bestSphere;
				bool collidedWithPolygon = true;

				// try ray - primitive intersection for all objects in the scene
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
						collidedWithPolygon = false;
					}
				}

				// check if any intersection has been found
				if (bestInt.distance < INFINITY) {
					// change the orientation of the ray to save computational time
					// in phong model
					r.dir = r.dir * -1;

					if (collidedWithPolygon) {
						// ray collided with a polygon first
						if (!bestPolygon.matType) {
							// the material of the polygon is a default material (Material class)
							Vertex pixelColor = phong(r, bestInt, materials[bestPolygon.matIdx], -1);

							setPixel(x, y, pixelColor.m_data[0], pixelColor.m_data[1], pixelColor.m_data[2]);
						}
						else {
							// the material of the polygon is a emissive material (EmissiveMaterial class)
							Vertex pixelColor = phong(r, bestInt, materials[bestPolygon.matIdx],-1);// TODO - phong with emissive material

							setPixel(x, y, pixelColor.m_data[0], pixelColor.m_data[1], pixelColor.m_data[2]);
						}
					}
					else {
						//  ray collided with a sphere first

						Vertex pixelColor = phong(r, bestInt, materials[bestSphere.matIdx], bestSphere.center.m_data[0]);

						setPixel(x, y, pixelColor.m_data[0], pixelColor.m_data[1], pixelColor.m_data[2]);
					}					
				}
				else {
					// No intersection with the ray and scene - set background color
					setPixel(x, y, clearColor.r, clearColor.g, clearColor.b);
				}
			}
		}
	}

	/// Loads 3 vertices from the vertex buffer and stores them as polygons in the scene
	void storePolygons() {
		
		Polygon p(
			vertexBuffer[0],
			vertexBuffer[1],
			vertexBuffer[2]
		);

		// get normal of the polygon
		cross(minus(p.b, p.a), minus(p.c, p.a), p.normal);
		normalize(p.normal);
		p.matIdx = (addingEmissiveMaterial) ? emmisiveMaterials.size() - 1 : materials.size() - 1;
		p.matType = addingEmissiveMaterial;
		polygons.push_back(p);
	}
};