//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2016/10/24
// Author: Jaroslav Sloup
// Author: Jaroslav Krivanek, Jiri Bittner, CTU Prague
// Edited: Jakub Hendrich, Daniel Meister, CTU Prague
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES

#include "sgl.h"
#include "sglhelper.h"
#include "matrix.hpp"

#include <math.h>
#include <vector>



/// Current error code.
static sglEErrorCode _libStatus = SGL_NO_ERROR;

/// Other structures
struct colorPixel {
	float r;
	float g;
	float b;
};

struct Vertex {
	float x;
	float y;
	float z;
	float w;

	Vertex() {};
	Vertex(float x, float y, float z, float w) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
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
	std::vector<Vertex> vertexBuffer;

	Context(int width, int height, int index) {
		this->width = width;
		this->height = height;
		this->index = index;
		colorBuffer = new float[width * height * 3];
		viewport.width = 1;
		viewport.height = 4;
		depthBuffer = new float[width * height];
		//TODO malloc vertexBuffer
		depthTest = false;
		clearColor.r = 0;
		clearColor.g = 0;
		clearColor.b = 0;
		drawingColor.r = 1;
		drawingColor.g = 1;
		drawingColor.b = 1;
		viewport.makeIdentity();
	}

	~Context() {

		delete[] colorBuffer;
		delete[] depthBuffer;

		/*for (unsigned int j = 0; j < vertexBuffer.size(); ++j)
			delete vertexBuffer[j];*/


		vertexBuffer.shrink_to_fit();
		projectionMatricesStack.clear();
		modelViewMatricesStack.clear();
		projectionMatricesStack.shrink_to_fit();
		modelViewMatricesStack.shrink_to_fit();
	}
};

/// SGL variables
bool sglBeginEndRunning = false;
bool sglBeginEndSceneRunning = false;
int maxStackSize = 100;
std::vector<Context*> contexts;
int contextCounter = 0;
Context* currentContext;
bool popFlag = false;


static inline void setErrCode(sglEErrorCode c)
{
	if (_libStatus == SGL_NO_ERROR)
		_libStatus = c;
}

//---------------------------------------------------------------------------
// sglGetError()
//---------------------------------------------------------------------------
sglEErrorCode sglGetError(void)
{
	sglEErrorCode ret = _libStatus;
	_libStatus = SGL_NO_ERROR;
	return ret;
}

//---------------------------------------------------------------------------
// sglGetErrorString()
//---------------------------------------------------------------------------
const char* sglGetErrorString(sglEErrorCode error)
{
	static const char *errStrigTable[] =
	{
	  "Operation succeeded",
	  "Invalid argument(s) to a call",
	  "Invalid enumeration argument(s) to a call",
	  "Invalid call",
	  "Quota of internal resources exceeded",
	  "Internal library error",
	  "Matrix stack overflow",
	  "Matrix stack underflow",
	  "Insufficient memory to finish the requested operation"
	};

	if ((int)error < (int)SGL_NO_ERROR || (int)error >(int)SGL_OUT_OF_MEMORY)
		return "Invalid value passed to sglGetErrorString()";

	return errStrigTable[(int)error];
}

//---------------------------------------------------------------------------
// Initialization functions
//---------------------------------------------------------------------------

void sglInit(void) {
	currentContext = 0;
	testFunction();
	try {
		// TODO ...
		// Alloc all needed data structures
	}
	catch (std::bad_alloc &e) {
		std::cerr << "Allocation failed! Message: " << e.what() << std::endl;
		_libStatus = SGL_OUT_OF_MEMORY;
		return;
	}
}

void sglFinish(void) {
	for (unsigned int i = 0; i < contexts.size(); i++) {
		Context *curr = contexts[i];
		delete curr;
	}

}

int sglCreateContext(int width, int height) {

	if (contextCounter > 50) { _libStatus = SGL_OUT_OF_RESOURCES; return -1; }

	Context *thisContext;
	try {
		thisContext = new Context(width, height, contextCounter);
	}
	catch (std::bad_alloc &e) {
		std::cerr << "Allocation failed! Message: " << e.what() << std::endl;
		_libStatus = SGL_OUT_OF_MEMORY;
		return -1;
	}

	contexts.push_back(thisContext);
	contextCounter++;
	return thisContext->index;
}

void sglDestroyContext(int id) {

	if (currentContext->index == id) { _libStatus = SGL_INVALID_OPERATION; return; }

	bool found = false;
	for (unsigned int i = 0; i < contexts.size(); i++) {
		int cont = contexts[i]->index;
		if (cont == id) {
			delete contexts[i];
			//contexts.erase(contexts.begin() + id);
			--contextCounter;
			found = true;
			break;
		}
	}
	if (!found) _libStatus = SGL_INVALID_VALUE;
}

void sglSetContext(int id) {

	bool found = false;
	for (unsigned int i = 0; i < contexts.size(); i++) {
		if (id == contexts[i]->index) {
			currentContext = contexts[i];
			found = true;
			break;
		}
	}
	if (!found) _libStatus = SGL_INVALID_VALUE;
}

int sglGetContext(void) {
	if (!currentContext) { _libStatus = SGL_INVALID_OPERATION; return -1; }
	return currentContext->index;
}

float *sglGetColorBufferPointer(void) {
	return (float*)currentContext->colorBuffer;
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor(float r, float g, float b, float alpha) {
	if (!currentContext || sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->clearColor.r = r;
	currentContext->clearColor.g = g;
	currentContext->clearColor.b = b;
}

void sglClear(unsigned what) {
	if (!currentContext || sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }
	//printf("Yavolal se\n");
	if (what == SGL_COLOR_BUFFER_BIT) {
		for (int i = 0; i < currentContext->width * currentContext->height; i += 3) {
			currentContext->colorBuffer[i] = currentContext->clearColor.r;
			currentContext->colorBuffer[i + 1] = currentContext->clearColor.g;
			currentContext->colorBuffer[i + 2] = currentContext->clearColor.b;
		}
	}
	else if (what == SGL_DEPTH_BUFFER_BIT) {
		for (int i = 0; i < currentContext->width * currentContext->height; i += 1) {
			currentContext->depthBuffer[i] = 1000000;
		}
	}
	else {
		_libStatus = SGL_INVALID_VALUE;
		return;
	}
}

void sglBegin(sglEElementType mode) {
	if (sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (mode > 8 || mode < 1) { _libStatus = SGL_INVALID_ENUM; return; }

	sglBeginEndRunning = true;
	currentContext->elementType = mode;
}

void setPixel(int x, int y) {
	if (y >= currentContext->height || y < 0 || x < 0 || x > currentContext->width) {
		return;
	}

	int position = x + y * currentContext->width;
	position *= 3;
	currentContext->colorBuffer[position] = currentContext->drawingColor.r;
	currentContext->colorBuffer[position + 1] = currentContext->drawingColor.g;
	currentContext->colorBuffer[position + 2] = currentContext->drawingColor.b;

}

void setSymetricalPixels(int x, int y, int xs, int ys)
{
	setPixel(xs + x, ys + y);
	setPixel(xs + y, ys + x);
	setPixel(xs + y, ys - x);
	setPixel(xs + x, ys - y);
	setPixel(xs - x, ys - y);
	setPixel(xs - y, ys - x);
	setPixel(xs - y, ys + x);
	setPixel(xs - x, ys + y);
}

void bresenhamCircle(float xs, float ys, float zs, float r) {

	if (currentContext->areaMode == SGL_POINT) {
		sglBegin(SGL_POINTS);
		sglVertex3f(xs, ys, zs);
		sglEnd();
	}

	Matrix matrix = currentContext->projectionMatricesStack.back() * currentContext->modelViewMatricesStack.back();

	Matrix stred(4, 1);
	float b[] = { xs, ys, zs, 1 };
	stred.initData(b);
	Matrix res = (matrix * stred) * (1 / stred.m_data[3][0]);
	int stx = currentContext->viewport.m_data[0][0] * res.m_data[0][0] + currentContext->viewport.m_data[2][0];
	int sty = currentContext->viewport.m_data[1][0] * res.m_data[1][0] + currentContext->viewport.m_data[3][0];

	Matrix MV = currentContext->modelViewMatricesStack.back();
	Matrix P = currentContext->projectionMatricesStack.back();

	float MVscale = MV.m_data[0][0] * MV.m_data[1][1] - MV.m_data[1][0] * MV.m_data[0][1];
	float Pscale = P.m_data[0][0] * P.m_data[1][1] - P.m_data[1][0] * P.m_data[0][1];
	r *= sqrt(MVscale * Pscale * currentContext->viewportScale);

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

void drawPoints()
{
	auto projSize = currentContext->projectionMatricesStack.size() - 1;
	auto mvSize = currentContext->modelViewMatricesStack.size() - 1;
	Matrix matrix = currentContext->projectionMatricesStack[projSize] * currentContext->modelViewMatricesStack[mvSize];

	for (Vertex vert : currentContext->vertexBuffer) {
		
		Matrix bod(4, 1);
		float b[] = { vert.x, vert.y, vert.z, vert.w };
		bod.initData(b);
		Matrix res = (matrix * bod) * (1 / bod.m_data[3][0]);

		int x = static_cast<int>(currentContext->viewport.m_data[0][0] * res.m_data[0][0] + currentContext->viewport.m_data[2][0]);
		int y = static_cast<int>(currentContext->viewport.m_data[1][0] * res.m_data[1][0] + currentContext->viewport.m_data[3][0]);
		
		int pointSize = currentContext->pointSize;
		
		
		if (pointSize % 2 == 0) {
			x = x - pointSize / 2 - 1;
			y = y - pointSize / 2 - 1;
		}
		else {
			x = x - (pointSize - 1) / 2;
			y = y - (pointSize - 1) / 2;
		}
		
		
		for (int a = 0; a < pointSize; a++) {
			for (int b = 0; b < pointSize; b++) {
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

void bresenhamLine(int x1, int x2, int y1, int y2)
{
	int c0, c1, p;

	if (x2 - x1 <= 0 && y2 - y1 <= 0) {
		int tempX = x1;
		int tempY = y1;

		x1 = x2;
		y1 = y2;

		x2 = tempX;
		y2 = tempY;

	}

	
	if(abs(y2 - y1) > abs(x2 - x1)) {//svisl�
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
	}else{
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

void approximationEllipse(float x, float y, float z, float a, float b) {

	if (currentContext->areaMode == SGL_POINT) {
		sglBegin(SGL_POINTS);
		sglVertex3f(x, y, z);
		sglEnd();
	}

	float x1 = a;
	float y1 = 0;
	float x2 = 0;
	float y2 = 0; 
	float alpha = 2 * M_PI / 40;
	float angle = 0;

	

	float CA = cos(alpha);
	float SA = sin(alpha);

	sglBegin(SGL_LINE_LOOP); //comment to use bresenham

	sglVertex3f(x + a * cos(angle), y + b * sin(angle), z);
	for (int i = 1; i < 40; i++) {
		//x2 = CA*x1 * a - SA*y1 * a;
		//y2 = SA*x1 * b  + CA*y1 * b;
		x2 = a * cos(i * alpha);
		y2 = b * sin(i * alpha);
		sglVertex3f(x + x2, y + y2, z);
		//bresenhamLine(x + x1, x + x2, y + y1, y + y2);

		//x1 = x2;
		//y1 = y2;

	}
	sglEnd();
}

void approximationArc(float x, float y, float z, float radius, float from, float to) {

	float x1;
	float y1;
	float x2;
	float y2;
	int alpha = 40.0 * abs(to - from) / (2 * M_PI);

	//float CA = cos(alpha);
	//float SA = sin(alpha);


	
	sglBegin(SGL_LINE_STRIP);

	

	//float num = 40 * abs(to - from) / (2 * M_PI);
	float step = abs(to - from) / alpha;

	
	x1 = radius * cos(from);
	y1 = radius * sin(from);

	sglVertex3f(x + x1, y + y1, z);

	for (int i = 0; i <= alpha; i++) {
		//x2 = CA * x1 - SA * y1;
		//y2 = SA * x1 + CA * y1;
		//bresenhamLine(x + x1, x + x2, y + y1, y + y2);
		//x1 = x2;
		//y1 = y2;

		//x2 = radius * cos(from + i * alpha);
		//y2 = radius * sin(from + i * alpha);
		x2 = radius * cos(from + i*step);
		y2 = radius * sin(from + i*step);
		//bresenhamLine(x + x1, x + x2, y + y1, y + y2);
		sglVertex3f(x + x2, y + y2, z);
		//x1 = x2;
		//y1 = y2;
		
	}
	sglEnd();
}


void drawLines()
{
	auto projSize = currentContext->projectionMatricesStack.size() - 1;
	auto mvSize = currentContext->modelViewMatricesStack.size() - 1;
	Matrix matrix = currentContext->projectionMatricesStack[projSize] * currentContext->modelViewMatricesStack[mvSize];
	//TODO postupn� projdeme v�echny vrcholy ve vertexBufferu a na ka�dou dal�� �se�ku zavol�me bresenhamLine();
	for (unsigned int i = 0; i < currentContext->vertexBuffer.size(); i += 2) {
		Vertex v1 = currentContext->vertexBuffer[i];
		Vertex v2 = currentContext->vertexBuffer[i + 1];
		Matrix bod1(4, 1);
		Matrix bod2(4, 1);
		float b1[] = { v1.x, v1.y, v1.z, v1.w };
		float b2[] = { v2.x, v2.y, v2.z, v2.w };
		bod1.initData(b1);
		bod2.initData(b2);
		Matrix res1 = (matrix * bod1) * (1 / bod1.m_data[3][0]);
		Matrix res2 = (matrix * bod2) * (1 / bod2.m_data[3][0]);

		int x1 = currentContext->viewport.m_data[0][0] * res1.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y1 = currentContext->viewport.m_data[1][0] * res1.m_data[1][0] + currentContext->viewport.m_data[3][0];
		int x2 = currentContext->viewport.m_data[0][0] * res2.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y2 = currentContext->viewport.m_data[1][0] * res2.m_data[1][0] + currentContext->viewport.m_data[3][0];

		/* int x1 = currentContext->vertexBuffer.at(i)->x;
		 int y1 = currentContext->vertexBuffer.at(i)->y;
		 int x2 = currentContext->vertexBuffer.at(i + 1)->x;
		 int y2 = currentContext->vertexBuffer.at(i + 1)->y;*/
		bresenhamLine(x1, x2, y1, y2);
	}
}

void drawLineStrip()
{
	auto projSize = currentContext->projectionMatricesStack.size() - 1;
	auto mvSize = currentContext->modelViewMatricesStack.size() - 1;
	Matrix matrix = currentContext->projectionMatricesStack[projSize] * currentContext->modelViewMatricesStack[mvSize];

	for (unsigned int i = 0; i < currentContext->vertexBuffer.size() - 1; i++) {

		Vertex v1 = currentContext->vertexBuffer[i];
		Vertex v2 = currentContext->vertexBuffer[i + 1];
		Matrix bod1(4, 1);
		Matrix bod2(4, 1);
		float b1[] = { v1.x, v1.y, v1.z, v1.w };
		float b2[] = { v2.x, v2.y, v2.z, v2.w };
		bod1.initData(b1);
		bod2.initData(b2);
		Matrix res1 = (matrix * bod1) * (1 / bod1.m_data[3][0]);
		Matrix res2 = (matrix * bod2) * (1 / bod2.m_data[3][0]);

		int x1 = currentContext->viewport.m_data[0][0] * res1.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y1 = currentContext->viewport.m_data[1][0] * res1.m_data[1][0] + currentContext->viewport.m_data[3][0];
		int x2 = currentContext->viewport.m_data[0][0] * res2.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y2 = currentContext->viewport.m_data[1][0] * res2.m_data[1][0] + currentContext->viewport.m_data[3][0];
		bresenhamLine(x1, x2, y1, y2);
	}

}

void drawLineLoop()
{

	auto projSize = currentContext->projectionMatricesStack.size() - 1;
	auto mvSize = currentContext->modelViewMatricesStack.size() - 1;
	Matrix matrix = currentContext->projectionMatricesStack[projSize] * currentContext->modelViewMatricesStack[mvSize];
	
	Vertex vert = currentContext->vertexBuffer[0];
	Matrix bod(4, 1);
	float b[] = { vert.x, vert.y, vert.z, vert.w };
	bod.initData(b);
	Matrix res = (matrix * bod) * (1 / bod.m_data[3][0]);

	int startx = currentContext->viewport.m_data[0][0] * res.m_data[0][0] + currentContext->viewport.m_data[2][0];
	int starty = currentContext->viewport.m_data[1][0] * res.m_data[1][0] + currentContext->viewport.m_data[3][0];

	int length = 0;

	for (unsigned int i = 0; i < currentContext->vertexBuffer.size() - 1; i++) {
		Vertex v1 = currentContext->vertexBuffer[i];
		Vertex v2 = currentContext->vertexBuffer[i + 1];
		Matrix bod1(4, 1);
		Matrix bod2(4, 1);
		float b1[] = { v1.x, v1.y, v1.z, v1.w };
		float b2[] = { v2.x, v2.y, v2.z, v2.w };
		bod1.initData(b1);
		bod2.initData(b2);
		Matrix res1 = (matrix * bod1) * (1 / bod1.m_data[3][0]);
		Matrix res2 = (matrix * bod2) * (1 / bod2.m_data[3][0]);

		int x1 = currentContext->viewport.m_data[0][0] * res1.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y1 = currentContext->viewport.m_data[1][0] * res1.m_data[1][0] + currentContext->viewport.m_data[3][0];
		int x2 = currentContext->viewport.m_data[0][0] * res2.m_data[0][0] + currentContext->viewport.m_data[2][0];
		int y2 = currentContext->viewport.m_data[1][0] * res2.m_data[1][0] + currentContext->viewport.m_data[3][0];
		bresenhamLine(x1, x2, y1, y2);
		length++;
	}

	Vertex vert2 = currentContext->vertexBuffer[length];
	Matrix bod2(4, 1);
	float b2[] = { vert2.x, vert2.y, vert2.z, vert2.w };
	bod2.initData(b2);
	Matrix res2 = (matrix * bod2) * (1 / bod2.m_data[3][0]);

	int endx = currentContext->viewport.m_data[0][0] * res2.m_data[0][0] + currentContext->viewport.m_data[2][0];
	int endy = currentContext->viewport.m_data[1][0] * res2.m_data[1][0] + currentContext->viewport.m_data[3][0];

	bresenhamLine(endx, startx, endy, starty);
}


void sglEnd(void) {
	if (!sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	switch (currentContext->elementType)
	{
	case SGL_POINTS: {
		drawPoints();
		//std::cout << "velikost bufferu" << currentContext->vertexBuffer.size() << std::endl;
		break;
	}
	case SGL_LINES: {
		drawLines();
		break;
	}
	case SGL_LINE_STRIP: {
		drawLineStrip();
		break;
	}
	case SGL_LINE_LOOP: {
		drawLineLoop();
		break;
	}
	default:
		break;
	}

	sglBeginEndRunning = false;
	currentContext->vertexBuffer.clear();
}


void sglVertex4f(float x, float y, float z, float w) {
	currentContext->vertexBuffer.push_back(Vertex(x, y, z, w));
}

void sglVertex3f(float x, float y, float z) {
	currentContext->vertexBuffer.push_back(Vertex(x, y, z, 1));
}

void sglVertex2f(float x, float y) {
	currentContext->vertexBuffer.push_back(Vertex(x, y, 0, 1));
}

void sglCircle(float x, float y, float z, float radius) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (radius < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	bresenhamCircle(x, y, z, radius);

}

void sglEllipse(float x, float y, float z, float a, float b) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (a < 0 || b < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	approximationEllipse(x, y, z, a, b);
}

void sglArc(float x, float y, float z, float radius, float from, float to) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (radius < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	approximationArc(x, y, z, radius, from, to);

}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode(sglEMatrixMode mode) {
	if (mode > 2 || mode < 0) { _libStatus = SGL_INVALID_ENUM; return; }

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->matrixMode = mode;
}

void sglPushMatrix(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (currentContext->matrixMode == SGL_PROJECTION) {

		if (maxStackSize == currentContext->projectionMatricesStack.size()) { _libStatus = SGL_STACK_OVERFLOW; return; }

		Matrix m = currentContext->projectionMatricesStack[currentContext->projectionMatricesStack.size() - 1];
		Matrix newM = m;
		currentContext->projectionMatricesStack.push_back(newM);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {

		if (maxStackSize == currentContext->modelViewMatricesStack.size()) { _libStatus = SGL_STACK_OVERFLOW; return; }

		Matrix m = currentContext->modelViewMatricesStack[currentContext->modelViewMatricesStack.size() - 1];
		Matrix newM = m;
		currentContext->modelViewMatricesStack.push_back(newM);
	}
}

void sglPopMatrix(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (currentContext->matrixMode == SGL_PROJECTION) {
		if (currentContext->projectionMatricesStack.size() <= 1 && !popFlag) { 
			_libStatus = SGL_STACK_UNDERFLOW; return; }

		currentContext->projectionMatricesStack.pop_back();
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		if (currentContext->modelViewMatricesStack.size() <= 1 && !popFlag) { 
			_libStatus = SGL_STACK_UNDERFLOW; return; }

		currentContext->modelViewMatricesStack.pop_back();
	}
}

void sglLoadIdentity(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (currentContext->matrixMode == SGL_PROJECTION) {
		bool emptyStack = currentContext->projectionMatricesStack.size() == 0;
		Matrix m;
		if (!emptyStack) {
			m = currentContext->projectionMatricesStack.back();
		}
		m.makeIdentity();
		Matrix newM = m;
		if (!emptyStack) {
			popFlag = true;
			sglPopMatrix();
			popFlag = false;
		}
		/*std::cout << "sglLoadIdentity: projection \n";
		newM.print();*/
		currentContext->projectionMatricesStack.push_back(newM);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		bool emptyStack = currentContext->modelViewMatricesStack.size() == 0;
		Matrix m;
		if (!emptyStack)
			m = currentContext->modelViewMatricesStack.back();
		m.makeIdentity();
		Matrix newM = m;
		if (!emptyStack) {
			popFlag = true;
			sglPopMatrix();
			popFlag = false;
		}
		/*std::cout << "sglLoadIdentity: modelview \n";
		newM.print();*/
		currentContext->modelViewMatricesStack.push_back(newM);
	}
}

void sglLoadMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	Matrix toLoad(4, 4);
	toLoad.initData(matrix);
	popFlag = true;
	sglPopMatrix();
	popFlag = false;
	if (_libStatus == SGL_NO_ERROR) {
		if (currentContext->matrixMode == SGL_PROJECTION) {
			currentContext->projectionMatricesStack.push_back(toLoad);
		}
		else if (currentContext->matrixMode == SGL_MODELVIEW) {
			currentContext->modelViewMatricesStack.push_back(toLoad);
		}
	}
}

void sglMultMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	Matrix toMult(4, 4);
	toMult.initData(matrix);

	if (currentContext->matrixMode == SGL_PROJECTION) {
		Matrix m = currentContext->projectionMatricesStack[currentContext->projectionMatricesStack.size() - 1];
		Matrix newM = m * toMult;
		popFlag = true;
		sglPopMatrix();
		popFlag = false;
		/*std::cout << "sglMult: projection \n";
		newM.print();*/
		currentContext->projectionMatricesStack.push_back(newM);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		Matrix m = currentContext->modelViewMatricesStack[currentContext->modelViewMatricesStack.size() - 1];
		Matrix newM = m * toMult;
		popFlag = true;
		sglPopMatrix();
		popFlag = false;
		/*std::cout << "sglMult: modelview \n";
		newM.print();*/
		currentContext->modelViewMatricesStack.push_back(newM);
	}
}

void sglTranslate(float x, float y, float z) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float m[] = {
	1, 0, 0, x, // col 1
	0, 1, 0, y, // col 2
	0, 0, 1, z, // col 3
	0, 0, 0, 1 };// col 4
	sglMultMatrix(m);
}

void sglScale(float scalex, float scaley, float scalez) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float m[] = {
	scalex, 0, 0, 0, // col 1
	0, scaley, 0, 0, // col 2
	0, 0, scalez, 0, // col 3
	0, 0, 0, 1 };// col 4
	sglMultMatrix(m);
}

void sglRotate2D(float angle, float centerx, float centery) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float sinus = sin(angle);
	float cosinus = cos(angle);

	/*float t_mat[] = {
	1, 0, 0, centerx, // col 1
	0, 1, 0, centery, // col 2
	0, 0, 1, 0, // col 3
	0, 0, 0, 1 };// col 4

	float minus_t_mat[] = {
	1, 0, 0, -centerx, // col 1
	0, 1, 0, -centery, // col 2
	0, 0, 1, 0, // col 3
	0, 0, 0, 1 };// col 4

	float rotate[] = {
	cosinus, -sinus, 0, 0, // col 1
	sinus, cosinus, 0, 0, // col 2
	0, 0, 1, 0, // col 3
	0, 0, 0, 1 };// col 4*/

	float all[] = {
		cosinus, -sinus, 0, -centerx * cosinus + centery * sinus + centerx, // col 1
		sinus, cosinus, 0,  -centerx * sinus - centery * cosinus + centery, // col 2
		0, 0, 1, 0, // col 3
		0, 0, 0, 1 };// col 4

	//sglMultMatrix(minus_t_mat);
	//sglMultMatrix(rotate);
	//sglMultMatrix(t_mat);
	sglMultMatrix(all);
}

void sglRotateY(float angle) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float sinus = sin(angle);
	float cosinus = cos(angle);

	float rotate[] = {
	cosinus,	0,	-sinus, 0,		// col 1
	0,			1,		 0, 0,		// col 2
	sinus,		0, cosinus, 0,		// col 3
	0,			0,		 0, 1 };	// col 4

	sglMultMatrix(rotate);
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far) { _libStatus = SGL_INVALID_VALUE; return; }

	float ortho[] = {
	2 / (right - left), 0, 0, -(right + left) / (right - left),
	0, 2 / (top - bottom), 0, -(top + bottom) / (top - bottom),
	0, 0, -2 / (far - near),  -(far + near) / (far - near),
	0, 0, 0, 1
	};
	sglMultMatrix(ortho);
}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far || near < 0 || far < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	float flustrum[] = {
	2 * near / (right - left), 0, (right + left) / (right - left), 0,
	0, 2 * near / (top - bottom), (top + bottom) / (top - bottom), 0,
	0, 0, (far + near) / (far - near), 2 * far * near / (far - near),
	0, 0, -1, 0
	};
	sglMultMatrix(flustrum);

}

void sglViewport(int x, int y, int width, int height) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (width < 0 || height < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->viewport.m_data[0][0] = width / 2.0f;
	currentContext->viewport.m_data[1][0] = height / 2.0f;
	currentContext->viewport.m_data[2][0] = x + width / 2.0f;
	currentContext->viewport.m_data[3][0] = y + height / 2.0f;

	currentContext->viewportScale = (width*height) / 4;
}

//---------------------------------------------------------------------------
// Attribute functions
//---------------------------------------------------------------------------

void sglColor3f(float r, float g, float b) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->drawingColor.r = r;
	currentContext->drawingColor.g = g;
	currentContext->drawingColor.b = b;
}

void sglAreaMode(sglEAreaMode mode) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (mode > 2 || mode < 0) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->areaMode = mode;
}

void sglPointSize(float size) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (size < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->pointSize = size;
}

void sglEnable(sglEEnableFlags cap) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (cap != 1) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->depthTest = cap;
}

void sglDisable(sglEEnableFlags cap) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (cap != 1) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->depthTest = cap;
}

//---------------------------------------------------------------------------
// RayTracing oriented functions
//---------------------------------------------------------------------------

void sglBeginScene() {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	sglBeginEndSceneRunning = true;

}

void sglEndScene() {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	sglBeginEndSceneRunning = false;

}

void sglSphere(const float x,
	const float y,
	const float z,
	const float radius)
{
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglMaterial(const float r,
	const float g,
	const float b,
	const float kd,
	const float ks,
	const float shine,
	const float T,
	const float ior)
{

}

void sglPointLight(const float x,
	const float y,
	const float z,
	const float r,
	const float g,
	const float b)
{
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglRayTraceScene() {
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglRasterizeScene() {
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglEnvironmentMap(const int width,
	const int height,
	float *texels)
{
	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglEmissiveMaterial(const float r,
	const float g,
	const float b,
	const float c0,
	const float c1,
	const float c2)
{
	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}