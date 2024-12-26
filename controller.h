#ifndef _CONTROLLER_H 
#define _CONTROLLER_H

#include "core.h"
#include "object.h"
#include "arcball.h"

class Controller {
public:
	Controller(int argc,char **argv, const char *windowName);
	~Controller();

	void BeginLoop();
	void RegisterObject(Object* object);

	void Reset();
	void Render();
	void InitCamera();

	int GetActiveObject(int mx, int my);

	void Quit();

	// Event handlers
	void Resize(GLFWwindow *window, int x, int y);
	void Keyboard(GLFWwindow * window, int key, int scancode,int action, int mods);
	void MouseButton(GLFWwindow *window, int btn, int action, int mods);
	void MouseMotion(GLFWwindow *window, double x, double y);
	void MouseScroll(GLFWwindow *window, double x, double y);

private:
	// Window management
	std::string _windowName;
	std::stringstream _titleInfo;
	GLFWwindow *_windowHandle;
	int _winX, _winY;

	// Input
	bool _isLeftKeyPressed, _isCtrlPressed, _isRightKeyPressed, _isMiddleKeyPressed;
	bool _isLightOn;
	
	double _prevCursorX, _prevCursorY;

	Object** _objects;
	int _objectNo;
	int _maxObjectNo;
	int _activeObj;

	void _ComputeFPS();
	double _startTime;
	double _stopTime;
	double _elapsedTime;
	int _numOfFrame;
	double _fps;

	// Perspective controls
	GLfloat _FOV;		// Field of View Angle
	GLfloat _aspect;	// Aspect Ratio
	GLfloat _nearClip;	// Near clipping plane distance
	GLfloat _farClip;	// Far clipping plane distance
};

#endif
