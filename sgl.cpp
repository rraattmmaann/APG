//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2016/10/24
// Author: Jaroslav Sloup
// Author: Jaroslav Krivanek, Jiri Bittner, CTU Prague
// Edited: Jakub Hendrich, Daniel Meister, CTU Prague
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES

#define MAX_RT_RECURSION_DEPTH 8
#define THRESHOLD 0.01f
#define SAMPLES	16

#include "sgl.h"
#include "context.h"
#include <vector>
#include <math.h>

/// Current error code.
static sglEErrorCode _libStatus = SGL_NO_ERROR;

/// SGL variables

bool sglBeginEndRunning = false;
bool sglBeginEndSceneRunning = false;
bool popFlag = false;
unsigned int maxStackSize = 100;
// poznamka: zatim jsem vektor contextu nechal dynamicky alokovany
std::vector<Context*> contexts;
int contextCounter = 0;
Context* currentContext;

using iterContext = std::vector<Context*>::iterator;

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
	currentContext = NULL;
	try {
		// TODO ...
		// Nothing to alloc yet, maybe in later homework parts?
	}
	catch (std::bad_alloc &e) {
		std::cerr << "Allocation failed! Message: " << e.what() << std::endl;
		_libStatus = SGL_OUT_OF_MEMORY;
		return;
	}
}

void sglFinish(void) {
	for (unsigned int i = 0; i < contexts.size(); i++) 
		delete contexts[i];	
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
	contexts.emplace_back(thisContext);
	contextCounter++;
	return thisContext->index;
}

void sglDestroyContext(int id) {

	if (!currentContext) { _libStatus = SGL_INVALID_OPERATION; return; }

	bool found = false;

	for (unsigned int i = 0; i < contexts.size(); ++i) {
		if (contexts[i]->index == id) {
			delete contexts[i];
			--contextCounter;
			found = true;
			break;
		}
	}
	if (!found) _libStatus = SGL_INVALID_VALUE;
}

void sglSetContext(int id) {

	bool found = false;

	for (unsigned int i = 0; i < contexts.size(); ++i) {
		if (contexts[i]->index == id) {
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
	if (sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->clearColor.r = r;
	currentContext->clearColor.g = g;
	currentContext->clearColor.b = b;
}

void sglClear(unsigned what) {
	if (sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	_libStatus = currentContext->clearBuffer(what);
}

void sglBegin(sglEElementType mode) {
	if (sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (mode > SGL_AREA_LIGHT || mode < SGL_POINTS) { _libStatus = SGL_INVALID_ENUM; return; }

	sglBeginEndRunning = true;
	currentContext->elementType = mode;
}

void sglEnd(void) {
	if (!sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (sglBeginEndSceneRunning) {
		currentContext->storePolygons();
		sglBeginEndRunning = false;
		currentContext->vertexBuffer.clear();
		return;
	}

	switch (currentContext->elementType)
	{
	case SGL_POINTS: {
		currentContext->drawPoints();
		break;
	}
	case SGL_LINES: {
		currentContext->drawLines();
		break;
	}
	case SGL_LINE_STRIP: {
		currentContext->drawLineStrip();
		break;
	}
	case SGL_LINE_LOOP: {
		currentContext->drawLineLoop();
		break;
	}
	case SGL_POLYGON: {
		currentContext->drawPolygon();
		break;
	}
	case SGL_AREA_LIGHT: {
		currentContext->drawAreaLight();
		break;
	}
	default:
		break;
	}

	sglBeginEndRunning = false;
	currentContext->vertexBuffer.clear();
}

void sglVertex4f(float x, float y, float z, float w) {
	currentContext->vertexBuffer.emplace_back(Vertex(x, y, z, w));
}

void sglVertex3f(float x, float y, float z) {
	currentContext->vertexBuffer.emplace_back(Vertex(x, y, z, 1.0f));
}

void sglVertex2f(float x, float y) {
	currentContext->vertexBuffer.emplace_back(Vertex(x, y, 0.0f, 1.0f));
}

void sglCircle(float x, float y, float z, float radius) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (radius < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->bresenhamCircle(x, y, z, radius);
}

void sglEllipse(float x, float y, float z, float a, float b) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (a < 0 || b < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->approximationEllipse(x, y, z, a, b);
}

void sglArc(float x, float y, float z, float radius, float from, float to) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (to >= 2.0f * M_PI) {
		to = fmod(to, 2.0f * float(M_PI));

		if (to == 0.0f) to = 2.0f * float(M_PI);
	}

	if (from >= 2.0f * M_PI) from = fmod(from, 2.0f * float(M_PI));

	if (radius < 0 || to < from) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->approximationArc(x, y, z, radius, from, to);
}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode(sglEMatrixMode mode) {
	if (mode != SGL_MODELVIEW && mode != SGL_PROJECTION) { _libStatus = SGL_INVALID_ENUM; return; }

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->matrixMode = mode;
}

void sglPushMatrix(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (currentContext->matrixMode == SGL_PROJECTION) {

		if (maxStackSize == currentContext->projectionMatricesStack.size()) { _libStatus = SGL_STACK_OVERFLOW; return; }

		Matrix m = currentContext->projectionMatricesStack.back();
		currentContext->projectionMatricesStack.emplace_back(m);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {

		if (maxStackSize == currentContext->modelViewMatricesStack.size()) { _libStatus = SGL_STACK_OVERFLOW; return; }

		Matrix m = currentContext->modelViewMatricesStack.back();
		currentContext->modelViewMatricesStack.emplace_back(m);
	}
}

void sglPopMatrix(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (currentContext->matrixMode == SGL_PROJECTION) {
		if (currentContext->projectionMatricesStack.size() <= 1 && !popFlag) { _libStatus = SGL_STACK_UNDERFLOW; return; }

		currentContext->projectionMatricesStack.pop_back();
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		if (currentContext->modelViewMatricesStack.size() <= 1 && !popFlag) { _libStatus = SGL_STACK_UNDERFLOW; return; }

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
			m.makeIdentity();
		}
		if (!emptyStack) {
			popFlag = true;
			sglPopMatrix();
			popFlag = false;
		}
		currentContext->projectionMatricesStack.emplace_back(m);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		bool emptyStack = currentContext->modelViewMatricesStack.size() == 0;
		Matrix m;
		if (!emptyStack) {
			m = currentContext->modelViewMatricesStack.back();
			m.makeIdentity();
		}		
		if (!emptyStack) {
			popFlag = true;
			sglPopMatrix();
			popFlag = false;
		}
		currentContext->modelViewMatricesStack.emplace_back(m);
	}
}

void sglLoadMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	Matrix toLoad(4, 4);
	toLoad.initData(matrix);
	popFlag = true;
	sglPopMatrix();
	popFlag = false;
		
	if (currentContext->matrixMode == SGL_PROJECTION) {
		currentContext->projectionMatricesStack.emplace_back(toLoad);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		currentContext->modelViewMatricesStack.emplace_back(toLoad);
	}
}

void sglMultMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	Matrix toMult(4, 4);
	toMult.initData(matrix);

	if (currentContext->matrixMode == SGL_PROJECTION) {
		Matrix m = currentContext->projectionMatricesStack.back();
		popFlag = true;
		sglPopMatrix();
		popFlag = false;
		currentContext->projectionMatricesStack.emplace_back(m * toMult);
	}
	else if (currentContext->matrixMode == SGL_MODELVIEW) {
		Matrix m = currentContext->modelViewMatricesStack.back();
		popFlag = true;
		sglPopMatrix();
		popFlag = false;
		currentContext->modelViewMatricesStack.emplace_back(m * toMult);
	}
}

void sglTranslate(float x, float y, float z) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float m[] = {
	1, 0, 0, 0,		// col 1
	0, 1, 0, 0,		// col 2
	0, 0, 1, 0,		// col 3
	x, y, z, 1 };	// col 4
	sglMultMatrix(m);
}

void sglScale(float scalex, float scaley, float scalez) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float m[] = {
	scalex, 0, 0, 0,	// col 1
	0, scaley, 0, 0,	// col 2
	0, 0, scalez, 0,	// col 3
	0, 0, 0, 1 };		// col 4
	sglMultMatrix(m);
}

void sglRotate2D(float angle, float centerx, float centery) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float sinus = sin(angle);
	float cosinus = cos(angle);

	// Already multiplied transformation: x_res = T^-1 * R * T * x

	float all[] = {
		cosinus, sinus, 0, 0,	// col 1
		-sinus, cosinus, 0,  0, // col 2
		0, 0, 1, 0,				// col 3
		-centerx * cosinus + centery * sinus + centerx, -centerx * sinus - centery * cosinus + centery, 0, 1 // col 4
	};

	sglMultMatrix(all);
}

void sglRotateY(float angle) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	float sinus = sin(-angle);
	float cosinus = cos(-angle);

	float rotate[] = {
	cosinus,	0,	-sinus,  0,		// col 1
	0,			1,		 0, 0,		// col 2
	sinus,		0, cosinus, 0,		// col 3
	0,			0,		 0, 1 };	// col 4

	sglMultMatrix(rotate);
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far) { _libStatus = SGL_INVALID_VALUE; return; }

	float ortho[] = {
	2 / (right - left), 0, 0, 0,	// col 1
	0, 2 / (top - bottom), 0, 0,	// col 2
	0, 0, -2 / (far - near),  0,		// col 3
	-(right + left) / (right - left), -(top + bottom) / (top - bottom), -(far + near) / (far - near), 1 };	// col 4

	sglMultMatrix(ortho);
}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far || near < 0 || far < 0) { _libStatus = SGL_INVALID_VALUE; return; }	

	float flustrum[] = {
	2 * near / (right - left),		0,							0,	0,	// col 1
	0,								2 * near / (top - bottom),	0,	0,	// col 2
	(right + left) / (right - left),(top + bottom) / (top - bottom),	-(far + near) / (far - near),	-1,	// col 3
	0,								0,							-2 * far * near / (far - near),		0 }; // col 4

	currentContext->viewport.m_data[0][1] = (far + near) / (near - far);
	currentContext->viewport.m_data[1][1] = 2 * far * near / (near - far);

	sglMultMatrix(flustrum);
}

void sglViewport(int x, int y, int width, int height) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (width < 0 || height < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->setViewport(x, y, width, height);
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

	if (mode > SGL_FILL || mode < SGL_POINT) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->areaMode = mode;
}

void sglPointSize(float size) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (size < 0) { _libStatus = SGL_INVALID_VALUE; return; }

	currentContext->pointSize = size;
}

void sglEnable(sglEEnableFlags cap) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (cap != SGL_DEPTH_TEST) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->depthTest = true;
}

void sglDisable(sglEEnableFlags cap) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (cap != SGL_DEPTH_TEST) { _libStatus = SGL_INVALID_ENUM; return; }

	currentContext->depthTest = false;
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
	if (sglBeginEndRunning || !sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }


	currentContext->spheres.push_back(
		Sphere(
			Vertex(x, y, z, 1.0f),
			radius,
			currentContext->materials.size() - 1
		)
	);
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
	currentContext->materials.push_back(
		Material(r, g, b, kd, ks, shine, T, ior)
	);
	currentContext->addingEmissiveMaterial = false;
}

void sglPointLight(const float x,
	const float y,
	const float z,
	const float r,
	const float g,
	const float b)
{
	if (sglBeginEndRunning || !sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->lights.push_back(
		Light(
			Vertex(x, y, z, 1.0f),
			r,g,b
		)
	);
}

void sglRayTraceScene() {
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	currentContext->startRt();
}

void sglRasterizeScene() {
	if (sglBeginEndRunning || sglBeginEndSceneRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
}

void sglEnvironmentMap(const int width,
	const int height,
	float *texels)
{
	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
	
	currentContext->enviromentMaps.push_back(
		EnviromentMap(width, height, texels)
	);
}

void sglEmissiveMaterial(const float r,
	const float g,
	const float b,
	const float c0,
	const float c1,
	const float c2)
{
	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return;	}

	currentContext->emissiveMaterials.push_back(
		EmissiveMaterial(r,g,b,c0,c1,c2)
	);
	currentContext->addingEmissiveMaterial = true;
}