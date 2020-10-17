#pragma once

#include <iostream>
#include <vector>

#include "sgl.h"
#include "matrix.hpp"

using iter = std::vector<Matrix>::iterator;

struct colorPixel {
	float r;
	float g;
	float b;
};

/// Holds all variables need fot the current context
class Context {
public:
	/// WIdth and height of the canvas
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
	std::vector<Matrix> vertexBuffer;

	Context() {}
	
	Context(int width, int height, int index) {
		this->width = width;
		this->height = height;
		this->index = index;
		colorBuffer = new float[width * height * 3];
		viewport.width = 1;
		viewport.height = 4;
		viewport.makeIdentity();
		depthBuffer = new float[width * height];
		depthTest = false;
		clearColor.r = 0;
		clearColor.g = 0;
		clearColor.b = 0;
		drawingColor.r = 1;
		drawingColor.g = 1;
		drawingColor.b = 1;
		vertexBuffer.reserve(width*height);
		modelViewMatricesStack.reserve(32);
		projectionMatricesStack.reserve(10);
	}

	~Context() {
		if (colorBuffer) delete[] colorBuffer;
		if (depthBuffer) delete[] depthBuffer;

		vertexBuffer.shrink_to_fit();
		projectionMatricesStack.clear();
		modelViewMatricesStack.clear();
		projectionMatricesStack.shrink_to_fit();
		modelViewMatricesStack.shrink_to_fit();
	}

	/// Creates and adds a vertex with specified coordidinates to the vertex buffer
	///		@param x[in] vertex x coordinate
	///		@param y[in] vertex y coordinate
	///		@param z[in] vertex z coordinate
	///		@param w[in] vertex w coordinate
	void addVertex(float x, float y, float z, float w);

	/// Computes vertex positions for an Arc
	///		@param x [in] circle center x coordinate
	///		@param y [in] circle center y coordinate
	///		@param z [in] circle center z coordinate
	///		@param a[in] x semiaxis length
	///		@param b[in] y semiaxis length
	void approximationArc(float x, float y, float z, float radius, float from, float to);

	/// Computes vertex positions for an Ellipse
	///		@param x[in] ellipse center x coordinate
	///		@param y[in] ellipse center y coordinate
	///		@param z[in] ellipse center z coordinate(depth)
	///		@param a[in] x semiaxis length
	///		@param b[in] y semiaxis length
	void approximationEllipse(float x, float y, float z, float a, float b);

	/// Renders a circle with given radius on defined location
	///		@param xs[in] circle center x coordinate
	///		@param ys[in] circle center y coordinate
	///		@param zs[in] circle center z coordinate
	///		@param r[in] cirlce radius
	void bresenhamCircle(float xs, float ys, float zs, float r);

	/// Connects 2 given vertices with a line using the Bersenham line algorithm
	///		@param x1[in] vertex No. 1 x coordinate
	///		@param x2[in] vertex No. 2 x coordinate
	///		@param y1[in] vertex No. 1 y coordinate
	///		@param y2[in] vertex No. 2 y coordinate
	void bresenhamLine(int x1, int x2, int y1, int y2);

	/// Takes all vertices from the vertex buffer and handles their rendering
	void drawPoints();

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	void drawLines();

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	void drawLineStrip();

	/// Draws lines between each 2 vertices from the vertex buffer and handles their rendering
	/// and connects the first and the last vertex
	void drawLineLoop();

	/// Places currently used color to the specifed position in the color buffer
	///		@param x[in] pixel x coordinate
	///		@param y[in] pixel y coordinate
	void setPixel(int x, int y);

	/// Places currently used color to 8 specifed positions in the color buffer
	/// while rendering a circle
	///		@param x[in] current x axis offset
	///		@param y[in] current y axis offset
	///		@param xs[in] circle center x coordinate
	///		@param ys[in] circle center y coordinate
	void setSymetricalPixels(int x, int y, int xs, int ys);
	
	/// Sets current viewport matrix
	///		@param x[in] viewport origin x coordinate(with respect to the canvas)
	///		@param y[in] viewport origin y coordinate(with respect to the canvas)
	///		@param width[in] viewport width in pixels
	///		@param height[in] viewport height in pixels
	void setViewport(float x, float y, float width, float height);

	/// Clear the buffer specified by the input value, returns error otherwise
	///		@param what[in] type of buffer to clear
	///		@return SGL_NO_ERROR on successful clear, SGL_INVALID_VALUE otherwise
	sglEErrorCode clearBuffer(unsigned what);
};