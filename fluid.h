#ifndef _FLUID_H
#define _FLUID_H

#include "core.h"
#include "arcball.h"

#define DT  0.1				// time step
#define RES 40			    // box resolution
#define N ((RES)-2)			// valid simulation area
#define SIZE ((RES)*(RES)*(RES))

#undef _I
#define _I(x,y,z) (((x)*(RES)*(RES))+((y)*(RES))+(z))	//FIXME

#undef FOR_ALL_CELL
#define FOR_ALL_CELL for (int i=1; i<=(N); i++) {\
	for (int j=1; j<=(N); j++) {\
		for (int k=1; k<=(N); k++) {

#undef END_FOR
#define END_FOR }}}

#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)

#ifndef ALMOST_EQUAL
#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)
#endif

#ifndef _I
#define _I(x,y,z) (((x)*(_RES)*(_RES))+((y)*(_RES))+(z))	//FIXME
#endif

#ifndef FOR_ALL_CELL
#define FOR_ALL_CELL for (int i=1; i<=(_N); i++) {\
	for (int j=1; j<=(_N); j++) {\
		for (int k=1; k<=(_N); k++) {
#define END_FOR }}}
#endif

#define SLICE_NUM			64.0f

class Fluid
{
private:
	// texture data
	unsigned char* _textureData;
	// texture handle
	unsigned int _hTexture;				
	//lighting infomations
	Eigen::Vector3f _lightDir;
	int _rayTemplate[4096][3];
	float *_volumeData;

	GLfloat _cubeVertices[8][3];
	GLfloat _cubeEdges[12][2][3];

	// draw the slices. mvMatrix must be the MODELVIEW_MATRIX
	void DrawSlices(GLdouble mvMatrix[16]);

	// intersect a plane with the cube, helper function for DrawSlices()
	// plane equation is Ax + By + Cz + D = 0
	std::vector<Eigen::Vector3f> IntersectEdges(float A, float B, float C, float D);

	void GenerateRayTemplate(int edgeLen);
	void CastLight(int edgelen, const float* dens, unsigned char* intensity);
	inline void LightRay(int x, int y, int z, int n, float decay, 
			const float* dens, unsigned char* intensity);

	void InitGL();

	int _SIZE;		//size of volume data
	int _N;			//
	int _RES;		//
protected:
	float _buffers[10][SIZE];
	float *_density, *_densityTmp;			// density
	float *_velX, *_velXTmp;			// velocity in x direction
	float *_velY, *_velYTmp;			// velocity in y direction
	float *_velZ, *_velZTmp;			// velocity in z direction
	float _dt;


	//for rendering
	Eigen::Vector3f _lightPos;
	bool _isLightSelected;
	bool _isRendering;
	bool _isDrawSliceOutline;

	GLfloat _rotX;
	GLfloat _rotY;
	GLfloat _depth;

	//for event handling
	bool _isLeftKeyPressed, _isCtrlPressed, _isRightKeyPressed, _isMiddleKeyPressed;

	Arcball _arcball;

	GLFWwindow* _windowHandle;
	int _winX, _winY;
protected:
	// simulation methods
	// beware: i changed stam's implementation from a destiation, source ordering
	// of the parameters to source, destiation
	void AddSource(float* src, float *dst);
	void AddBuoyancy();
	void EnforceBoundary(int b, float* quantity);
	void Diffuse(int b, float* velTmp, float* vel, float diff);
	void Advect(int b, float* quantityTmp, float* quantity, float* velX, float* velY, float* velZ);
	void Project();
	void VorticityConfinement();

	void VelocityStep();
	void DensityStep();

	// utility methods
	void ClearBuffer(float* buf);
	void ClearSources(void);

	void GenerateSmoke();
	bool LightSelected(double mouseX, double mouseY);

public:
	float sd[SIZE], su[SIZE], sv[SIZE], sw[SIZE], sT[SIZE];	// sources for density and velocities
	float diffusion, viscosity, buoyancy, vc_eps;

	Fluid(void);
	~Fluid(void);

	virtual void SimulateStep();
	virtual void Show();
	const float* GetDensity();

	virtual void MouseButton(GLFWwindow *window, int button,int action,int mods);
	virtual void MouseMotion(GLFWwindow *window, double nx, double ny);
	virtual void Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods);

	virtual void Reset();
	virtual void Resize(GLFWwindow *window, int x, int y);
	virtual void MouseScroll(GLFWwindow *window, double nx, double ny);
	void RegisterParentWindow(GLFWwindow* windowHandle);

	void SetLightPostion(Eigen::Vector3f &pos);
	void SetRendering(bool isRendering);
	void SetSliceOutline(bool isDrawSliceOutline);
	void FillTexture();		// generate texture from smoke density 
	void Render();					// draw the volume
	// draw the outline of the cube
	void DrawCube();
	void DrawLight();
	void DrawVolumeData();
};

#endif

