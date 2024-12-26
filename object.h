#ifndef DRAWER_H 
#define DRAWER_H

#include "core.h"
#include "arcball.h"

////////////////////////////////////////////////////////////////////////////////

class Object: public EventListener{
public:
	Object();
	~Object();

	virtual void SimulateStep() {};
	virtual void Show();
	virtual void Reset();

	virtual void MouseButton(GLFWwindow *window, int button,int action,int mods);
	virtual void MouseMotion(GLFWwindow *window, double nx, double ny);
	virtual void MouseScroll(GLFWwindow *window, double nx, double ny);
	virtual void Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods) {}
	virtual void Resize(GLFWwindow *window, int x, int y);

	void RegisterParentWindow(GLFWwindow* windowHandle);
protected:
	std::fstream _file;

	GLfloat _rotX;
	GLfloat _rotY;
	GLfloat _depth;

	//for event handling
	bool _isLeftKeyPressed, _isCtrlPressed, _isRightKeyPressed, _isMiddleKeyPressed;

	Arcball _arcball;

	GLFWwindow* _windowHandle;
	int _winX, _winY;
	
	//for rendering
	float _stepSize;
};

////////////////////////////////////////////////////////////////////////////////

#endif
