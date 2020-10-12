//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2016/10/24
// Author: Jaroslav Sloup
// Author: Jaroslav Krivanek, Jiri Bittner, CTU Prague
// Edited: Jakub Hendrich, Daniel Meister, CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "sglhelper.h"

#include <vector>

/// Current error code.
static sglEErrorCode _libStatus = SGL_NO_ERROR;

/// Other structures
struct colorPixel {
	float r;
	float g;
	float b;
};

struct vertex {
    float x;
    float y;
    float z;
    float w;

    vertex() {};
    vertex(float x, float y) {
        this->x = x;
        this->y = y;
    }
    vertex(float x, float y, float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    vertex(float x, float y, float z, float w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }
};

struct context {
    int width, height;
    int index;
    colorPixel clearColor;

    float* colorBuffer;
    
    float* depthBuffer;
    float vieportMatrix[4][4];

    bool depthTest;
    enum sglEMatrixMode matrixMode;
    enum sglEElementType elementType;

    colorPixel drawingColor;
    float pointSize;

    std::vector<std::array<float, 16>> modelViewMatricesStack;
    std::vector<std::array<float, 16>> projectionMatricesStack;
    std::vector<vertex*> vertexBuffer;

    context(int width, int height, int index) {
        this->width = width;
        this->height = height;
        this->index = index;
        colorBuffer = (float*)malloc(width * height * 3 * sizeof(float));
        depthBuffer = (float*)malloc(width * height * sizeof(float));
        //TODO malloc vertexBuffer
        depthTest = false;
        clearColor.r = 0;
        clearColor.g = 0;
        clearColor.b = 0;
    }
};

/// SGL variables

bool sglBeginEndRunning = false;

std::vector<context*> contexts;
int contextCounter = 0;
context* currentContext;
std::array<float, 16> identity{ {1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0} };


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

  if ((int)error < (int)SGL_NO_ERROR || (int)error > (int)SGL_OUT_OF_MEMORY)
    return "Invalid value passed to sglGetErrorString()"; 

  return errStrigTable[(int)error];
}

//---------------------------------------------------------------------------
// Initialization functions
//---------------------------------------------------------------------------

void sglInit(void) {
	currentContext = 0;
	testFunction();
}

void sglFinish(void) {
    for (int i = 0; i < contexts.size(); i++) {
        delete contexts[i];
    }

}

int sglCreateContext(int width, int height) {
    context *thisContext = new context(width, height, contextCounter);
    contexts.push_back(thisContext);
    contextCounter++;
	return thisContext->index;
}

void sglDestroyContext(int id) {
    for (int i = 0; i < contexts.size(); i++) {
        int cont = contexts[i]->index;
        if (cont == id) {
            delete contexts[i];
            contexts.erase(contexts.begin() + id);
        }
    }
}

void sglSetContext(int id) {
    for (int i = 0; i < contexts.size(); i++) {
        if (id == contexts[i]->index) {
            currentContext = contexts[i];
            break;
        }
    }
}

int sglGetContext(void) {
    return currentContext->index;
}

float *sglGetColorBufferPointer(void) {
    return currentContext->colorBuffer;
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor(float r, float g, float b, float alpha) {
    std::cout << currentContext->index << std::endl;
    currentContext->clearColor.r = r;
    currentContext->clearColor.g = g;
    currentContext->clearColor.b = b;
}

void sglClear(unsigned what) {
    if (what == SGL_COLOR_BUFFER_BIT) {
        for (int i = 0; i < currentContext->width * currentContext->height; i += 3) {
            currentContext->colorBuffer[i] = currentContext->clearColor.r;
            currentContext->colorBuffer[i + 1] = currentContext->clearColor.g;
            currentContext->colorBuffer[i + 2] = currentContext->clearColor.b;
        }
    }
    else {
        for (int i = 0; i < currentContext->width * currentContext->height; i += 1) {
            currentContext->depthBuffer[i] = INFINITY;
        }
    }
}

void sglBegin(sglEElementType mode) {
    sglBeginEndRunning = true;
    currentContext->elementType = mode;
}

void sglEnd(void) {
    sglBeginEndRunning = false;
}

void sglVertex4f(float x, float y, float z, float w) {

}

void sglVertex3f(float x, float y, float z) {

}

void sglVertex2f(float x, float y) {
    currentContext->vertexBuffer.push_back(new vertex(x,y));
}

void sglCircle(float x, float y, float z, float radius) {

}

void sglEllipse(float x, float y, float z, float a, float b) {

}

void sglArc(float x, float y, float z, float radius, float from, float to) {

}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode(sglEMatrixMode mode) {
    currentContext->matrixMode = mode;
}

void sglPushMatrix(void) {
 
}

void sglPopMatrix(void) {

}

void sglLoadIdentity(void) {
    if (currentContext->matrixMode == SGL_PROJECTION) {
        currentContext->projectionMatricesStack.push_back(identity);
    }
    if (currentContext->matrixMode == SGL_MODELVIEW) {
        currentContext->modelViewMatricesStack.push_back(identity);
    }
}

void sglLoadMatrix(const float *matrix) {

}

void sglMultMatrix(const float *matrix) {

}

void sglTranslate(float x, float y, float z) {

}

void sglScale(float scalex, float scaley, float scalez) {

}

void sglRotate2D(float angle, float centerx, float centery) {

}

void sglRotateY(float angle) {

}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {

}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {

}

void sglViewport(int x, int y, int width, int height) {
    currentContext->vieportMatrix[0][0] = width/2;
    currentContext->vieportMatrix[0][1] = 0;
    currentContext->vieportMatrix[0][2] = 0;
    currentContext->vieportMatrix[0][3] = x+width/2;
    currentContext->vieportMatrix[1][0] = 0;
    currentContext->vieportMatrix[1][1] = height/2;
    currentContext->vieportMatrix[1][2] = 0;
    currentContext->vieportMatrix[1][3] = y + height / 2;
    currentContext->vieportMatrix[2][0] = 0;
    currentContext->vieportMatrix[2][1] = 0;
    currentContext->vieportMatrix[2][2] = 1;
    currentContext->vieportMatrix[2][3] = 0;
    currentContext->vieportMatrix[3][0] = 0;
    currentContext->vieportMatrix[3][1] = 0;
    currentContext->vieportMatrix[3][2] = 0;
    currentContext->vieportMatrix[3][3] = 1;
}

//---------------------------------------------------------------------------
// Attribute functions
//---------------------------------------------------------------------------

void sglColor3f(float r, float g, float b) {
    currentContext->drawingColor.r = r;
    currentContext->drawingColor.g = g;
    currentContext->drawingColor.b = b;

}

void sglAreaMode(sglEAreaMode mode) {

}

void sglPointSize(float size) {
    currentContext->pointSize = size;
}

void sglEnable(sglEEnableFlags cap) {
    currentContext->depthTest = cap;
}

void sglDisable(sglEEnableFlags cap) {
    currentContext->depthTest = cap;
}

//---------------------------------------------------------------------------
// RayTracing oriented functions
//---------------------------------------------------------------------------

void sglBeginScene() {

}

void sglEndScene() {

}

void sglSphere(const float x,
               const float y,
               const float z,
               const float radius) 
{

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

}

void sglRayTraceScene() {

}

void sglRasterizeScene() {

}

void sglEnvironmentMap(const int width,
                       const int height,
                       float *texels)
{

}

void sglEmissiveMaterial(const float r,
                         const float g,
                         const float b,
                         const float c0,
                         const float c1,
                         const float c2)
{

}