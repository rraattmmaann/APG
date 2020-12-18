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

	/// Vector containing all polygon in the scene
	std::vector<Polygon> polygons;

	/// Vector containing all point light sources in the scene
	std::vector<Light> lights;

	/// Vector containing all area light sources in the scene
	std::vector<Polygon> areaLights;

	/// Vector containing all emissive materials present in the scene
	std::vector<EmissiveMaterial> emissiveMaterials;

	/// Vector containing all environment maps present in the scene
	std::vector<EnviromentMap> enviromentMaps;

	/// Vector containing all spheres present in the scene
	std::vector<Sphere> spheres;

	/// Vector containing all materials present in the scene
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
		int steps = int(40.0f * (to - from) / (2.0f * M_PI));

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
		float alpha = 2.0f * float(M_PI) / 40.0f;
		float x1 = 1.0f;
		float y1 = 0.0f;
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
		stred = stred * (1.0f / stred.m_data[3]);
		int stx = int(viewport.m_data[0][0] * stred.m_data[0] + viewport.m_data[2][0]);
		int sty = int(viewport.m_data[1][0] * stred.m_data[1] + viewport.m_data[3][0]);

		// Compute the transformed radius
		float MVscale = MV.m_data[0][0] * MV.m_data[1][1] - MV.m_data[1][0] * MV.m_data[0][1];
		float Pscale = P.m_data[0][0] * P.m_data[1][1] - P.m_data[1][0] * P.m_data[0][1];
		r *= sqrt(MVscale * Pscale * viewportScale);
		
		if (areaMode == SGL_FILL) {

			int x1;
			int x2;

			int to = int(sty + round(r));
			float sqrtContent = (r*r);
			float outerContent = stx + 0.5f;

			for (int count = int(sty - r); count <= to; ++count) {
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
			y = int(r);
			p = int(3.0f - 2.0f * r);
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
			res = res * (1.0f / res.m_data[3]);

			int x = static_cast<int>(viewport.m_data[0][0] * res.m_data[0] + viewport.m_data[2][0]);
			int y = static_cast<int>(viewport.m_data[1][0] * res.m_data[1] + viewport.m_data[3][0]);

			int currPointSize = int(pointSize);

			if (currPointSize % 2 == 0) {
				x = int(x - currPointSize * 0.5f - 1.f);
				y = int(y - currPointSize * 0.5f - 1.f);
			}
			else {
				x = int(x - (currPointSize - 1.f) * 0.5f);
				y = int(y - (currPointSize - 1.f) * 0.5f);
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
			v1 = v1 * (1.0f / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1.0f / v2.m_data[3]);

			bresenhamLine(
				int(viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0]),
				int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0])
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
			v1 = v1 * (1.0f / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1.0f / v2.m_data[3]);

			bresenhamLine(
				int(viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0]),
				int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0])
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
		v = v * (1.0f / v.m_data[3]);

		int startx = int(viewport.m_data[0][0] * v.m_data[0] + viewport.m_data[2][0]);
		int starty = int(viewport.m_data[1][0] * v.m_data[1] + viewport.m_data[3][0]);

		int length = 0;

		for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {
			Vertex v1 = vertexBuffer[i];
			Vertex v2 = vertexBuffer[i + 1];
			v1 = PVM * v1;
			v1 = v1 * (1.0f / v1.m_data[3]);
			v2 = PVM * v2;
			v2 = v2 * (1.0f / v2.m_data[3]);

			bresenhamLine(
				int(viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0]),
				int(viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0]),
				int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0])
			);
			length++;
		}

		Vertex v2 = vertexBuffer[length];
		v2 = PVM * v2;
		v2 = v2 * (1.0f / v2.m_data[3]);

		int endx = int(viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0]);
		int endy = int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]);

		bresenhamLine(endx, startx, endy, starty);
	}

	/// Finds intersection in the x axis in height of y on line between vertices a and b
	///		@param y[in] current height of the scanline
	///		@param a[in] bottom end of the edge
	///		@param b[in] upper end of the edge
	///		@param intersections[in] vector of x-coordinates of all found intersections
	inline void getIntersections(int y, Vertex &a, Vertex &b, std::vector<int> &intersections) {
		int x1 = int(a.m_data[0]);
		int y1 = int(a.m_data[1]);
		int x2 = int(b.m_data[0]);		
		int y2 = int(b.m_data[1]);

		if (x2-x1 == 0) {
			intersections.emplace_back(x1);
			return;
		}

		float slope = (float)(y2 - y1) / (float)(x2 - x1);
		
		float shift = y1 - slope * x1;

		// y = slope * x + shift
		// dosadime za y promenou y xd
		int x = int((y - shift) / slope);
		intersections.emplace_back(x);
	}

	/// Finds intersection in the x axis in height of y on line between vertices a and b
	/// Computes the z value of the intersection as well
	///		@param y[in] current height of the scanline
	///		@param a[in] bottom end of the edge
	///		@param b[in] upper end of the edge
	///		@param intersections[in] vector of couple (x-coordinate, depth) of all found intersections
	inline void getIntersectionsDepth(int y, Vertex &a, Vertex &b, std::vector<std::tuple<int, float>> &intersections) {
		int x1 = int(a.m_data[0]);
		int y1 = int(a.m_data[1]);
		int x2 = int(b.m_data[0]);
		int y2 = int(b.m_data[1]);
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
		int x = int((y - shift) / slope);
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
				v1 = v1 * (1.0f / v1.m_data[3]);
				v2 = PVM * v2;
				v2 = v2 * (1.0f / v2.m_data[3]);

				int y1 = int(viewport.m_data[1][0] * v1.m_data[1] + viewport.m_data[3][0]);
				int y2 = int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]);
				
				screenSpaceVertices.emplace_back(Vertex(
					viewport.m_data[0][0] * v1.m_data[0] + viewport.m_data[2][0],
					float(y1),
					viewport.m_data[0][1] * v1.m_data[2] + viewport.m_data[1][1],
					1)
				);
				screenSpaceVertices.emplace_back(Vertex(
					viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0],
					float(y2), 
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
			v = v * (1.0f / v.m_data[3]);
			Vertex v2 = vertexBuffer[length];
			v2 = PVM * v2;
			v2 = v2 * (1.0f / v2.m_data[3]);

			int startx = int(viewport.m_data[0][0] * v.m_data[0] + viewport.m_data[2][0]);
			int starty = int(viewport.m_data[1][0] * v.m_data[1] + viewport.m_data[3][0]);
			float startz = viewport.m_data[0][1] * v.m_data[2] + viewport.m_data[1][1];		

			int endx = int(viewport.m_data[0][0] * v2.m_data[0] + viewport.m_data[2][0]);
			int endy = int(viewport.m_data[1][0] * v2.m_data[1] + viewport.m_data[3][0]);
			float endz = viewport.m_data[0][1] * v2.m_data[2] + viewport.m_data[1][1];

			screenSpaceVertices.emplace_back(Vertex(float(endx), float(endy), endz, 1));
			screenSpaceVertices.emplace_back(Vertex(float(startx), float(starty), startz, 1));

			// Run the scanline algorithm on the screenspace vertices
			for (int y = yMax; y > yMin; --y) {
				
				// Find all intersections of the current scanline with the edges
				// POZNAMKA: toto lze zrychlit predpocitanim hran a tedy i jejich smernic(kazda hrana je dana prinkou s predpisem y = ax + b)
				//           protoze vypocet se stale opakuje pro kazdou uroven scanline
				for (unsigned int i = 0; i < screenSpaceVertices.size(); i+=2) {

					Vertex a = screenSpaceVertices[i];
					Vertex b = screenSpaceVertices[i + 1];

					int y1 = int(a.m_data[1]);
					int y2 = int(b.m_data[1]);

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
								float z = 1.0f /((z2*j - z1*j + temp) / (x2 - x1));
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
	void setViewport(int x, int y, int width, int height) {
		viewport.m_data[0][0] = width / 2.0f;
		viewport.m_data[1][0] = height / 2.0f;
		viewport.m_data[2][0] = x + width / 2.0f;
		viewport.m_data[3][0] = y + height / 2.0f;

		viewportScale = (width*height) / 4.0f;
	}

	/// Clear the buffer specified by the input value, returns error otherwise
	///		@param what[in] type of buffer to clear
	///		@return SGL_NO_ERROR on successful clear, SGL_INVALID_VALUE otherwise
	sglEErrorCode clearBuffer(unsigned what) {
		if (what == (SGL_COLOR_BUFFER_BIT | SGL_DEPTH_BUFFER_BIT)) {
			for (int i = 0; i < width * height; ++i) {
				depthBuffer[i] = INFINITY;
			}
			for (int i = 0; i < width * height * 3; i += 3) {
				colorBuffer[i] = clearColor.r;
				colorBuffer[i + 1] = clearColor.g;
				colorBuffer[i + 2] = clearColor.b;
			}
		} else if (what == SGL_COLOR_BUFFER_BIT) {
			for (int i = 0; i < width * height * 3; i += 3) {
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
	///		@param ret[in] vector cointaining the resulting color of the pixel (r,g,b,1)
	void phong(Ray &ray, Intersection &intersection, Material &mat, Vertex &ret) {
		
		Vertex matColor(mat.r* mat.kd, mat.g* mat.kd, mat.b* mat.kd, 1);

		// Go through each light in the scene and sum their influence on given point
		for (unsigned int i = 0; i < lights.size(); i++) {
			Light l = lights[i];
			
			Vertex L = l.position - intersection.position;
			float dist = vectorLength(l.position, intersection.position); // distance from light to ray intersection point
			normalize(L);

			Ray r(intersection.position, L);
			
			// Check if the light source is shadowed by an object
			// if it is hidden, continue with next light source
			bool lightSourceHidden = false;

			for (unsigned int j = 0; j < polygons.size(); ++j) {

				Intersection Int = polygons[j].intersects(r);
				
				if (abs(dist - Int.distance) < 0.01f) continue; // eliminate floating point arithmetic inprecision
				
				if (Int.distance != INFINITY && Int.distance < dist) {
					lightSourceHidden = true;
					break;
				}
			}
			if (lightSourceHidden) continue;
				
			for (unsigned int j = 0; j < spheres.size(); ++j) {

				Intersection Int = spheres[j].intersects(r);
				
				if (abs(dist - Int.distance) < 0.01f) continue; // eliminate floating point arithmetic inprecision
				
				if (Int.distance != INFINITY && Int.distance < dist) {
					lightSourceHidden = true;
					break;
				}
			}
			if (lightSourceHidden)	continue;

			
			Vertex light(l.r, l.g, l.b, 1.0f);			
			Vertex minusL = L * -1.0f;
			Vertex R = minusL - intersection.normal * dot(minusL, intersection.normal) * 2.0f;

			Vertex v = matColor * light * std::max(0.0f, dot(intersection.normal, L));			// diffuse
			v += light * std::pow(std::max(0.0f, dot(R, ray.dir)), mat.shine) * mat.ks;			// specular

			ret += v;
		}
	}

	/// Phong illumination model modified to work with area lights in the scene
	///		@param ray[in] ray casted from the viewer(camera) to the primitive
	///		@param intersection[in] the point of intersection of the ray with the primitive
	///		@param mat[in] material used for the lightning model
	///		@param ret[in] vector cointaining the resulting color of the pixel (r,g,b,1)
	void areaPhong(Ray &ray, Intersection &intersection, Material &mat, Vertex &ret) {

		// Color of material of the object hit by given ray
		Vertex matColor(mat.r* mat.kd, mat.g* mat.kd, mat.b* mat.kd, 1);		
		Vertex p, L, zero{ 0,0,0,0 };
		float r1, r2, A, u, v, w, dist;
		bool lightSourceHidden;


		// Go through each area light in the scene and sum their influence on given point
		for (unsigned int i = 0; i < areaLights.size(); i++) {

			Vertex e1 = areaLights[i].b - areaLights[i].a;		// 1st side of the triangle
			Vertex e2 = areaLights[i].c - areaLights[i].a;		// 2nd side of the triangle
			cross(e1, e2, p);
			A = vectorLength(zero, p) * 0.5f;
			EmissiveMaterial eM = emissiveMaterials[areaLights[i].matIdx];	// emissive material of the triangle
			Vertex lightColor(eM.r, eM.g, eM.b, 1.0f);							// Color of the light (emisssive triangle)

			// Cast 16 sample rays from given hitpoint
			for (int j = 0; j < SAMPLES; ++j) {

				// 2 random numbers in <0,1>
				r1 = ((float)rand() / (RAND_MAX));
				r2 = ((float)rand() / (RAND_MAX));

				// Compute numbers in approopriate interval
				
				if (r1 + r2 > 1) {
					u = 1.0f - r1;
					v = 1.0f - r2;
				}
				else {
					u = r1;
					v = r2;
				}

				// Create random point light in the emissive triangle 
				p = areaLights[i].a + e1 * u + e2 * v;

				// Define vector to the created point light
				L = p - intersection.position;
				dist = vectorLength(p, intersection.position); // distance from light to ray intersection point
				normalize(L);

				// Create vector to the created point light
				Ray r(intersection.position, L);

				// Check if the light source is shadowed by an object
				// if it is hidden, continue with next light source
				lightSourceHidden = false;

				for (unsigned int j = 0; j < polygons.size(); ++j) {

					Intersection Int = polygons[j].intersects(r);

					if (abs(dist - Int.distance) < 0.01f) continue; // eliminate floating point arithmetic inprecision

					if (Int.distance != INFINITY && Int.distance < dist) {
						lightSourceHidden = true;
						break;
					}
				}
				if (lightSourceHidden) continue;

				Vertex minusL = L * -1.0f;

				// Weight of this point should influence the final color
				w = (A * dot(areaLights[i].normal, minusL)) / (SAMPLES * (eM.c0 + dist * eM.c1 + dist * dist * eM.c2));

				// diffuse only, multiply by weight
				ret += matColor * lightColor * std::max(0.0f, dot(intersection.normal, L)) * w; 
			}
		}
	}

	/// Accepts given ray, determines the color of point in the scene
	/// and return given color, alternatively sends reflected ray
	///		@param r[in] ray to be traced
	///		@param depth[in] current depth of the ray recursion within the scene
	///		@return Vertex containing a color of hit pixel
	Vertex traceRay(Ray &r, int depth) {

		Intersection	bestInt;
		Polygon			bestPolygon;
		Sphere			bestSphere;
		Material		bestMaterial;
		EmissiveMaterial bestEmMaterial;
		bool			collidedWithPolygon = true;

		// try ray - primitive intersection for all objects in the scene
		for (unsigned int i = 0; i < polygons.size(); ++i) {
			Intersection Int = polygons[i].intersects(r);
			if (Int.distance < bestInt.distance) {
				bestInt = Int;
				bestPolygon = polygons[i];
				if (!bestPolygon.matType)
					bestMaterial = materials[polygons[i].matIdx];
				else
					bestEmMaterial = emissiveMaterials[polygons[i].matIdx];
			}
		}

		for (unsigned int i = 0; i < spheres.size(); ++i) {
			Intersection Int = spheres[i].intersects(r);
			if (Int.distance < bestInt.distance) {
				bestInt = Int;
				bestSphere = spheres[i];
				bestMaterial = materials[spheres[i].matIdx];
				collidedWithPolygon = false;
			}
		}

		Vertex hitColor;
		
		// check if any intersection has been found
		if (bestInt.distance < INFINITY) {

			// change the orientation of the ray to save computational time in phong model
			Ray temp(r.origin, r.dir * -1.0f);

			if (collidedWithPolygon) {				// ray collided with a polygon first				
				if (!bestPolygon.matType) {			// the material of the polygon is a default material (Material class)
					if (areaLights.size() != 0)
						areaPhong(temp, bestInt, bestMaterial, hitColor);	// areaLights
					else
						phong(temp, bestInt, bestMaterial, hitColor);		// point lights
				}				
					
				else			
					// the material of the polygon is an emissive material (EmissiveMaterial class)
					// return the color of the emissive material light -> it was hit and nothing reflects
					return Vertex(bestEmMaterial.r, bestEmMaterial.g, bestEmMaterial.b,1);	
			}
			else {
				// ray collided with a sphere first
				phong(temp, bestInt, bestMaterial, hitColor);
			}
		}
		else {
			// No intersection with the ray and scene - set background color or env. map

			if (enviromentMaps.size() > 0) {
				
				// An enviroment map is defined, compute the pixel color from the texture
				float c = sqrt(std::pow(r.dir.m_data[0],2) + std::pow(r.dir.m_data[1], 2));
				float val = 0.0f;				
				if (c > 0.0f) {
					val = acos(r.dir.m_data[2])/(2 * c * (float)M_PI);
				}
				
				// Get normalized U and V coordinates in the texture
				// Note that 'V' has to be flipped -> the buffer is filled from TOP left corner,
				// whereas the origin of u,v coords is in (0,0) on BOTTOM left
				float u = 0.5f + val * r.dir.m_data[0];
				float v = 0.5f + val * r.dir.m_data[1];

				// Get real X,Y coordinates according to env. map width and height
				int x = int(u * enviromentMaps[0].width);
				int y = int((1-v) * enviromentMaps[0].height);

				// Compute the index of the pixel in the texel array
				int idx = (y * enviromentMaps[0].width + x) * 3;

				// Return RGB color
				return Vertex(
					enviromentMaps[0].texels[idx],
					enviromentMaps[0].texels[idx + 1],
					enviromentMaps[0].texels[idx + 2],
					1
				);
			
			} else
				// Return clearColor value - no enviroment map
				return Vertex(clearColor.r, clearColor.g, clearColor.b, 1);
		}

		if (bestMaterial.ks == 0.0f && bestMaterial.T == 0.0f)
			return hitColor;				// material specular coef is 0 -> the ray will not reflect IF the material is not transparent		

		// Compute reflected ray
		Ray reflectedRay(bestInt.position, (r.dir) - (bestInt.normal * 2 * dot(bestInt.normal, r.dir)));
		normalize(reflectedRay.dir);

		if (bestMaterial.T != 0.0f) { 
			// material is transparent - we need to generate refracted rays
			Ray refractedRay;
			
			Vertex normal = bestInt.normal;
			float gamma;
			float d = dot(r.dir, normal);
			if (d < 0.0f) {
				gamma = 1.0f / bestMaterial.ior;
				refractedRay.refracted = true;
			} else {
				gamma = bestMaterial.ior;
				d = -d;
				normal = normal * -1;
				refractedRay.refracted = false;
			}
			float sqrterm = 1.0f - gamma * gamma * (1.0f - d * d);
			if (sqrterm > 0.0f) {
				sqrterm = d * gamma + sqrt(sqrterm);
				
				refractedRay.origin = bestInt.position;
				refractedRay.dir = normal * -sqrterm + r.dir * gamma;				
				normalize(refractedRay.dir);

				hitColor += traceRay(refractedRay, depth + 1) * bestMaterial.T;
			}
		}
				
		return (depth + 1 > MAX_RT_RECURSION_DEPTH) ?
			// We reached max recursion depth - return current material color only 
			hitColor : 

			// Max depth has not been reached - combine pixel color with output of further recursion
			(hitColor + traceRay(reflectedRay, depth + 1) * bestMaterial.ks);
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
		// define the origin of all rays
		Ray r;
		r.origin = invMV * Vertex(0, 0, 0, 1);
		
		// go through all of the pixels and cast rays from the camera through each pixel
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				
				// direction of current ray				
				r.dir = MVP * Vertex(float(x) , float(y), -1.0f, 1.0f) - r.origin;
				normalize(r.dir);

				// Trace this ray and get the color of pixel on (x,y)
				Vertex pixelColor = traceRay(r, 0);
				setPixel(x, y, pixelColor.m_data[0], pixelColor.m_data[1], pixelColor.m_data[2]);
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
		cross(p.b - p.a, p.c - p.a, p.normal);
		normalize(p.normal);
		p.matIdx = (addingEmissiveMaterial) ? emissiveMaterials.size() - 1 : materials.size() - 1;
		p.matType = addingEmissiveMaterial;
		polygons.push_back(p);

		if (addingEmissiveMaterial)	
			areaLights.push_back(p);			
	}
};