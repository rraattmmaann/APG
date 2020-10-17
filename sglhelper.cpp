#define _USE_MATH_DEFINES
#include <cmath> 
#include <vector>

#include "sglhelper.h"

void Context::addVertex(float x, float y, float z, float w) {
	Matrix m(4, 1);
	float b[] = { x,y,z,w };
	m.initData(b, false);
	vertexBuffer.emplace_back(m);
}

void Context::approximationArc(float x, float y, float z, float radius, float from, float to) {

	float x1;
	float y1;
	float x2;
	float y2;
	int steps = 40 * (to - from) / (2 * M_PI);

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
		//bresenhamLine(x + x1, x + x2, y + y1, y + y2);
		//x1 = x2;
		//y1 = y2;
		//x2 = radius * cos(from + i * alpha);
		//y2 = radius * sin(from + i * alpha);
		//x2 = radius * cos(from + i*step);
		//y2 = radius * sin(from + i*step);
		//bresenhamLine(x + x1, x + x2, y + y1, y + y2);
		sglVertex3f(x + x2, y + y2, z);
		x1 = x2;
		y1 = y2;
	}
	sglEnd();
}

void Context::approximationEllipse(float x, float y, float z, float a, float b) {

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

void Context::bresenhamCircle(float xs, float ys, float zs, float r) {

	if (areaMode == SGL_POINT) {
		sglBegin(SGL_POINTS);
		sglVertex3f(xs, ys, zs);
		sglEnd();
	}
	Matrix MV = modelViewMatricesStack.back();
	Matrix P = projectionMatricesStack.back();
	Matrix matrix = P * MV;

	Matrix stred(4, 1);
	float b[] = { xs, ys, zs, 1 };
	stred.initData(b, true);
	Matrix res = (matrix * stred) * (1 / stred.m_data[3][0]);
	int stx = viewport.m_data[0][0] * res.m_data[0][0] + viewport.m_data[2][0];
	int sty = viewport.m_data[1][0] * res.m_data[1][0] + viewport.m_data[3][0];

	float MVscale = MV.m_data[0][0] * MV.m_data[1][1] - MV.m_data[1][0] * MV.m_data[0][1];
	float Pscale = P.m_data[0][0] * P.m_data[1][1] - P.m_data[1][0] * P.m_data[0][1];
	r *= sqrt(MVscale * Pscale * viewportScale);

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

void Context::bresenhamLine(int x1, int x2, int y1, int y2) {
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

void Context::drawPoints() {
	Matrix matrix = projectionMatricesStack.back() * modelViewMatricesStack.back();

	for (Matrix vert : vertexBuffer) {
		Matrix res = (matrix * vert) * (1 / vert.m_data[3][0]);


		int x = static_cast<int>(viewport.m_data[0][0] * res.m_data[0][0] + viewport.m_data[2][0]);
		int y = static_cast<int>(viewport.m_data[1][0] * res.m_data[1][0] + viewport.m_data[3][0]);

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


		//for (int a = static_cast<int>(x - pointSize / 2) ; a < static_cast<int>(x + pointSize / 2); a++) {
		//	for (int b = static_cast<int>(y - pointSize / 2) ; b < static_cast<int>(y + pointSize / 2); b++) {
		//		setPixel(a, b);
		//	}
		//}
	}
}

void Context::drawLines() {
	Matrix matrix = projectionMatricesStack.back() * modelViewMatricesStack.back();
	
	for (unsigned int i = 0; i < vertexBuffer.size(); i += 2) {
		Matrix v1 = vertexBuffer[i];
		Matrix v2 = vertexBuffer[i + 1];
		Matrix res1 = (matrix * v1) * (1 / v1.m_data[3][0]);
		Matrix res2 = (matrix * v2) * (1 / v2.m_data[3][0]);

		bresenhamLine(
			viewport.m_data[0][0] * res1.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[0][0] * res2.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[1][0] * res1.m_data[1][0] + viewport.m_data[3][0],
			viewport.m_data[1][0] * res2.m_data[1][0] + viewport.m_data[3][0]
		);
	}
}

void Context::drawLineStrip() {
	Matrix matrix = projectionMatricesStack.back() * modelViewMatricesStack.back();

	for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {

		Matrix v1 = vertexBuffer[i];
		Matrix v2 = vertexBuffer[i + 1];
		Matrix res1 = (matrix * v1) * (1 / v1.m_data[3][0]);
		Matrix res2 = (matrix * v2) * (1 / v2.m_data[3][0]);

		bresenhamLine(
			viewport.m_data[0][0] * res1.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[0][0] * res2.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[1][0] * res1.m_data[1][0] + viewport.m_data[3][0],
			viewport.m_data[1][0] * res2.m_data[1][0] + viewport.m_data[3][0]
		);
	}

}

void Context::drawLineLoop() {
	Matrix matrix = projectionMatricesStack.back() * modelViewMatricesStack.back();

	Matrix vert = vertexBuffer[0];
	Matrix res = (matrix * vert) * (1 / vert.m_data[3][0]);

	int startx = viewport.m_data[0][0] * res.m_data[0][0] + viewport.m_data[2][0];
	int starty = viewport.m_data[1][0] * res.m_data[1][0] + viewport.m_data[3][0];

	int length = 0;

	for (unsigned int i = 0; i < vertexBuffer.size() - 1; i++) {
		Matrix v1 = vertexBuffer[i];
		Matrix v2 = vertexBuffer[i + 1];
		Matrix res1 = (matrix * v1) * (1 / v1.m_data[3][0]);
		Matrix res2 = (matrix * v2) * (1 / v2.m_data[3][0]);

		bresenhamLine(
			viewport.m_data[0][0] * res1.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[0][0] * res2.m_data[0][0] + viewport.m_data[2][0],
			viewport.m_data[1][0] * res1.m_data[1][0] + viewport.m_data[3][0],
			viewport.m_data[1][0] * res2.m_data[1][0] + viewport.m_data[3][0]
		);
		length++;
	}

	Matrix vert2 = vertexBuffer[length];
	/*Matrix bod2(4, 1);
	float b2[] = { vert2.x, vert2.y, vert2.z, vert2.w };
	bod2.initData(b2);*/
	Matrix res2 = (matrix * vert2) * (1 / vert2.m_data[3][0]);

	int endx = viewport.m_data[0][0] * res2.m_data[0][0] + viewport.m_data[2][0];
	int endy = viewport.m_data[1][0] * res2.m_data[1][0] + viewport.m_data[3][0];

	bresenhamLine(endx, startx, endy, starty);
}

void Context::setPixel(int x, int y) {
	if (y >= height || y < 0 || x < 0 || x > width) {
		return;
	}

	int position = x + y * width;
	position *= 3;
	colorBuffer[position] = drawingColor.r;
	colorBuffer[position + 1] = drawingColor.g;
	colorBuffer[position + 2] = drawingColor.b;
}

void Context::setSymetricalPixels(int x, int y, int xs, int ys) {
	setPixel(xs + x, ys + y);
	setPixel(xs + y, ys + x);
	setPixel(xs + y, ys - x);
	setPixel(xs + x, ys - y);
	setPixel(xs - x, ys - y);
	setPixel(xs - y, ys - x);
	setPixel(xs - y, ys + x);
	setPixel(xs - x, ys + y);
}

void Context::setViewport(float x, float y, float width, float height) {
	viewport.m_data[0][0] = width / 2.0f;
	viewport.m_data[1][0] = height / 2.0f;
	viewport.m_data[2][0] = x + width / 2.0f;
	viewport.m_data[3][0] = y + height / 2.0f;

	viewportScale = (width*height) / 4;
}

sglEErrorCode Context::clearBuffer(unsigned what) {
	if (what == SGL_COLOR_BUFFER_BIT) {
		for (int i = 0; i < width * height; i += 3) {
			colorBuffer[i] = clearColor.r;
			colorBuffer[i + 1] = clearColor.g;
			colorBuffer[i + 2] = clearColor.b;
		}
	}
	else if (what == SGL_DEPTH_BUFFER_BIT) {
		for (int i = 0; i < width * height; i += 1) {
			depthBuffer[i] = 1000000;
		}
	}
	else {
		return SGL_INVALID_VALUE;
	}
	return SGL_NO_ERROR;
}