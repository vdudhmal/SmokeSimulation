#include "core.h"
#include "fluid.h"

#define WINDOW_NAME		 "Smoke3D"

//#define VERTEX_PROG_NAME1   "backface.vert"
//#define FRAGMENT_PROG_NAME1 "backface.frag"
//#define VERTEX_PROG_NAME2   "raycasting.vert"
//#define FRAGMENT_PROG_NAME2 "raycasting.frag"

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

void _ComputeFPS()
{
	_numOfFrame++;
	_stopTime = glfwGetTime();
	_elapsedTime += (_stopTime - _startTime);
	_startTime = glfwGetTime();
	if(_elapsedTime > 1) {
		_fps = _numOfFrame / _elapsedTime;
		//std::cout << _fps  << " "<< _numOfFrame << std::endl;
		_elapsedTime = 0;
		_numOfFrame = 0;
	}
	
	//clear stringstream buffer
	_titleInfo.str("");
	_titleInfo << _windowName;
	_titleInfo << "     FPS: ";
	_titleInfo << std::setprecision(4) << _fps;
	_titleInfo.width(2);
}

//set callback functions
static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void resize(GLFWwindow *window, int x,int y)
{
	_winX = x;
	_winY = y;
	
	for(int i = 0 ; i < _objectNo ; i++)
		_objects[i]->Resize(window, x, y);
}

static void keyboard(GLFWwindow * window, int key, int scancode,int action, int mods)		
{	
	if (action == GLFW_PRESS) {
		switch(key) {
			case GLFW_KEY_ESCAPE:		// Escape
				glFinish();
				glfwDestroyWindow(_windowHandle);
				exit(0);
				break;
			case GLFW_KEY_R:			//reset
				#ifdef DEBUG_LEVEL
					std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
				#endif
				glfwGetWindowSize(_windowHandle, &_winX, &_winY);
				_FOV = 45.0f;
				_aspect = (GLfloat)_winX / _winY;
				_nearClip = 0.1f;
				_farClip = 100.0f;

				#if 1
					glEnable(GL_DEPTH_TEST);
				#else
					glDisable(GL_DEPTH_TEST);
				#endif
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();

				// Set perspective projection
				gluPerspective(_FOV, _aspect, _nearClip, _farClip);

				glViewport(0, 0, _winX, _winY);

				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				for(int i = 0 ; i < _objectNo ; i++)
					_objects[i]->Reset();
				break;
			case GLFW_KEY_LEFT_CONTROL:
				_isCtrlPressed = true;
				break;
		}
	}
	else if(action == GLFW_RELEASE)
		switch (key) {
			case GLFW_KEY_LEFT_CONTROL:
				_isCtrlPressed = false;
				break;
		}

	//route message to active object
	_objects[_activeObj]->Keyboard(window, key, scancode, action, mods);
}

static void mousebutton(GLFWwindow *window,int button,int action,int mods)	
{
	//get active object and then transfer the message to the object
	if(action == GLFW_PRESS) {
		glfwGetCursorPos(window, &_prevCursorX, &_prevCursorY);
	}
	_objects[_activeObj]->MouseButton(window, button, action, mods);
}

static void mousemotion(GLFWwindow *window, double x, double y)
{
	_objects[_activeObj]->MouseMotion(window, x, y);
}

static void mousescroll(GLFWwindow *window, double x, double y)
{
	_objects[_activeObj]->MouseScroll(window, x, y);
}

int main(int argc, char **argv) {
	_windowName = WINDOW_NAME;

	//for timer
	_elapsedTime = 0;
	_numOfFrame = 0;
	_fps = 0;

	_winX = 1280;
	_winY = 960;

	_activeObj = -1;
	_objectNo = 0;
	_maxObjectNo = 100;
	_objects = new Object*[_maxObjectNo];

	_isCtrlPressed = _isLeftKeyPressed = _isMiddleKeyPressed = _isRightKeyPressed = false;
	_isLightOn = true;
	_prevCursorX = _prevCursorY = 0;

	// Initialize components
	if (!glfwInit()) {
		std::cerr << "glfwInit() failed!" << std::endl;
		exit(-1);
	}

	// Create the window
	_windowHandle = glfwCreateWindow(_winX, _winY, WINDOW_NAME, NULL, NULL);
	if (!_windowHandle) {
		std::cerr << "Create Window failed"  << std::endl;
		exit(-1);
	}
	glfwMakeContextCurrent(_windowHandle);
	glfwSetWindowPos(_windowHandle, 0, 0);

	// Background color
	glClearColor( 0., 0., 0., 1. );

	// Callbacks
	glfwSetErrorCallback(error_callback);
	glfwSetMouseButtonCallback(_windowHandle, mousebutton);
	glfwSetScrollCallback(_windowHandle, mousescroll);
	glfwSetCursorPosCallback(_windowHandle, mousemotion);
	glfwSetKeyCallback(_windowHandle, keyboard);
	glfwSetWindowSizeCallback(_windowHandle, resize);

	if(glewInit() != GLEW_OK) {
		std::cerr << "glewInit() failed!" << std::endl;
		exit(-1);
	}

	Object* object = new Fluid;

#if 0
	//enable shader
	object->RegisterShader(VERTEX_PROG_NAME1, GL_VERTEX_SHADER);
	object->RegisterShader(FRAGMENT_PROG_NAME1, GL_FRAGMENT_SHADER);
	object->RegisterShader(VERTEX_PROG_NAME2, GL_VERTEX_SHADER);
	object->RegisterShader(FRAGMENT_PROG_NAME2, GL_FRAGMENT_SHADER);
	object->EnableShader();
#endif

		#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif

	glfwGetWindowSize(_windowHandle, &_winX, &_winY);
	_FOV = 60.0f;
	_aspect = (GLfloat)_winX / _winY;
	_nearClip = 0.1f;
	_farClip = 100.0f;

#if 1
	glEnable(GL_DEPTH_TEST);
#else
	glDisable(GL_DEPTH_TEST);
#endif
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set perspective projection
	gluPerspective(_FOV, _aspect, _nearClip, _farClip);

	glViewport(0, 0, _winX, _winY);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	object->RegisterParentWindow(_windowHandle);
	object->Reset();
	_objects[_objectNo++] = object;
	_activeObj++;
	
	while (!glfwWindowShouldClose(_windowHandle))
    {
    	/* Render here */
		//compute fps and display 
		_ComputeFPS();
		glfwSetWindowTitle(_windowHandle, _titleInfo.str().c_str() );

		//Begin drawing scene
		glfwGetWindowSize(_windowHandle, &_winX, &_winY);
		_FOV = 45.0f;
		_aspect = (GLfloat)_winX / _winY;
		_nearClip = 0.1f;
		_farClip = 100.0f;

	#if 1
		glEnable(GL_DEPTH_TEST);
	#else
		glDisable(GL_DEPTH_TEST);
	#endif
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		// Set perspective projection
		gluPerspective(_FOV, _aspect, _nearClip, _farClip);

		glViewport(0, 0, _winX, _winY);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		for(int i = 0 ; i < _objectNo ; i++) {
			_objects[i]->SimulateStep();
			_objects[i]->Show();
		}

        /* Swap front and back buffers */
        glfwSwapBuffers(_windowHandle);

        /* Poll for and process events */
        glfwPollEvents();
    }

	glFinish();
	glfwDestroyWindow(_windowHandle);
	glfwTerminate();

	return 0;
}
