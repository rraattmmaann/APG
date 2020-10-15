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

struct Context {
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

    Context(int width, int height, int index) {
        this->width = width;
        this->height = height;
        this->index = index;
        colorBuffer = new float[width * height * 3];
        
        depthBuffer = new float[width * height];
        //TODO malloc vertexBuffer
        depthTest = false;
        clearColor.r = 0;
        clearColor.g = 0;
        clearColor.b = 0;
        drawingColor.r = 1;
        drawingColor.g = 1;
        drawingColor.b = 1;
    }
};

/// SGL variables

bool sglBeginEndRunning = false;
bool sglBeginEndSceneRunning = false;

std::vector<Context*> contexts;
int contextCounter = 0;
Context* currentContext;
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

		free(curr->colorBuffer);
		free(curr->depthBuffer);

		for (unsigned int j = 0; j < curr->vertexBuffer.size(); ++j) {
			delete curr->vertexBuffer[i];
		}

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
            contexts.erase(contexts.begin() + id);
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
    if (found) {
        //std::cout << id << std::endl;
    }
	if (!found) _libStatus = SGL_INVALID_VALUE;
}

int sglGetContext(void) {
	if (!currentContext) { _libStatus = SGL_INVALID_OPERATION; return -1; }
    return currentContext->index;
}

float *sglGetColorBufferPointer(void) {
    //std::cout << currentContext->index << std::endl;
    //std::cout << currentContext->colorBuffer[100 + 100 * currentContext->width] << " returning buffer" << std::endl;
    printf("color buffer: %f\n ", currentContext->colorBuffer[3*(100 + 100 * currentContext->width)]);
    return (float*)currentContext->colorBuffer;
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor(float r, float g, float b, float alpha) {
	if (!currentContext || sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }
   
	//std::cout << currentContext->index << std::endl;
    currentContext->clearColor.r = r;
    currentContext->clearColor.g = g;
    currentContext->clearColor.b = b;
}

void sglClear(unsigned what) {
	if (!currentContext || sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }
    printf("Yavolal se\n");
	if (what == SGL_COLOR_BUFFER_BIT) {
        for (int i = 0; i < currentContext->width * currentContext->height; i += 3) {
            currentContext->colorBuffer[i] = currentContext->clearColor.r;
            currentContext->colorBuffer[i + 1] = currentContext->clearColor.g;
            currentContext->colorBuffer[i + 2] = currentContext->clearColor.b;
        }
    }
    else if (what == SGL_DEPTH_BUFFER_BIT){
        for (int i = 0; i < currentContext->width * currentContext->height; i += 1) {
            currentContext->depthBuffer[i] = INFINITY;
        }
	}
	else {
		_libStatus = SGL_INVALID_VALUE;
		return;
	}
}

void sglBegin(sglEElementType mode) {
	if (sglBeginEndRunning) { _libStatus = SGL_INVALID_OPERATION; return; }
	
	if (mode > 8 || mode < 1) {_libStatus = SGL_INVALID_ENUM; return;}
    
	sglBeginEndRunning = true;
    currentContext->elementType = mode;
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
}

void drawPoints()
{

    for (vertex *vert : currentContext->vertexBuffer) {
        float x = currentContext->vieportMatrix[0][0] * vert->x + currentContext->vieportMatrix[0][2];
        float y = currentContext->vieportMatrix[0][1] * vert->y + currentContext->vieportMatrix[0][3];

        for (int a = x; a < x + currentContext->pointSize; a++) {
            for (int b = y; b < y + currentContext->pointSize; b++) {
                setPixel(a, b);
            }
        }
    
    }
    
}

void drawLines()
{
    //TODO postupn� projdeme v�echny vrcholy ve vertexBufferu a na ka�dou dal�� �se�ku zavol�me bresenhamLine();

}

void drawLineStrip()
{
    //TODO
}

void drawLineLoop()
{
    //TODO
}

void bresenhamLine(int x1, int x2, int y1, int y2)
{
    int c0, c1, p;
    c0 = 2 * (y2 - y1);
    c1 = c0 - 2 * (x2 - x1);
    p = c0 - (x2 - x1);
    for (int i = x1 + 1; i <= x2; i++) {
        if (p < 0) {
            p += c0;
        }
        else {
            p += c1;
            y1++;
        }
        setPixel(i, y1);
    }

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

void sglVertex4f(float x, float y, float z, float w) {

}

void sglVertex3f(float x, float y, float z) {

}

void sglVertex2f(float x, float y) {
    currentContext->vertexBuffer.push_back(new vertex(x,y));
}

void sglCircle(float x, float y, float z, float radius) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (radius < 0) { _libStatus = SGL_INVALID_VALUE; return; }

}

void sglEllipse(float x, float y, float z, float a, float b) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (a < 0 || b < 0) { _libStatus = SGL_INVALID_VALUE; return; }
}

void sglArc(float x, float y, float z, float radius, float from, float to) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (radius < 0) { _libStatus = SGL_INVALID_VALUE; return; }
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
 
}

void sglPopMatrix(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglLoadIdentity(void) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }
	
    if (currentContext->matrixMode == SGL_PROJECTION) {
        currentContext->projectionMatricesStack.push_back(identity);
    }
    if (currentContext->matrixMode == SGL_MODELVIEW) {
        currentContext->modelViewMatricesStack.push_back(identity);
    }
}

void sglLoadMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglMultMatrix(const float *matrix) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglTranslate(float x, float y, float z) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglScale(float scalex, float scaley, float scalez) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglRotate2D(float angle, float centerx, float centery) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglRotateY(float angle) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far) { _libStatus = SGL_INVALID_VALUE; return; }


}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (left == right || top == bottom || near == far || near < 0 || far < 0) { _libStatus = SGL_INVALID_VALUE; return; }

}

void sglViewport(int x, int y, int width, int height) {

	if (sglBeginEndRunning || contextCounter < 1) { _libStatus = SGL_INVALID_OPERATION; return; }

	if (width < 0 || height < 0) { _libStatus = SGL_INVALID_VALUE; return; }
    
    currentContext->vieportMatrix[0][0] = width/2.0f;
    currentContext->vieportMatrix[0][1] = height/2.0f;
    currentContext->vieportMatrix[0][2] = x + width / 2.0f;
    currentContext->vieportMatrix[0][3] = (y+height)/2.0f;
   
    
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