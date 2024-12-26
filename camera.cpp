#include "core.h"
#include "camera.h"

Camera::Camera(GLFWwindow* windowHandle) {
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	_windowHandle = windowHandle;

	Reset();
}

void Camera::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
}

void Camera::Reset()
{

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

}

void Camera::SetAspect(GLfloat aspect)		
{
	_aspect = aspect;
}
