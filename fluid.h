#ifndef _FLUID_H
#define _FLUID_H

#include "core.h"
#include "arcball.h"
#include "renderer.h"

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

class Renderer;

class Fluid
{
protected:
	float _buffers[10][SIZE];
	float *_density, *_densityTmp;			// density
	float *_velX, *_velXTmp;			// velocity in x direction
	float *_velY, *_velYTmp;			// velocity in y direction
	float *_velZ, *_velZTmp;			// velocity in z direction
	float _dt;


	//for rendering
	Renderer * _renderer;				//register a renderer
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
};

#endif

