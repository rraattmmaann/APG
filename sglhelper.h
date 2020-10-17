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

struct Context {
	int width, height;
	int index;
	colorPixel clearColor;

	float* colorBuffer;

	float* depthBuffer;
	Matrix viewport;
	float viewportScale;

	bool depthTest;
	enum sglEMatrixMode matrixMode;
	enum sglEElementType elementType;
	enum sglEAreaMode areaMode;

	colorPixel drawingColor;
	float pointSize;

	std::vector<Matrix> modelViewMatricesStack;
	std::vector<Matrix> projectionMatricesStack;
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

	void approximationEllipse(float x, float y, float z, float a, float b);
	void setPixel(int x, int y);
	void setSymetricalPixels(int x, int y, int xs, int ys);
	void bresenhamCircle(float xs, float ys, float zs, float r);
	void bresenhamLine(int x1, int x2, int y1, int y2);
	void drawPoints();
	void approximationArc(float x, float y, float z, float radius, float from, float to);
	void drawLines();
	void drawLineStrip();
	void drawLineLoop();
	void addVertex(float x, float y, float z, float w);
	void setViewport(float x, float y, float width, float height);
	sglEErrorCode clearBuffer(unsigned what);
};

