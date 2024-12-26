#include "core.h"
#include "object.h"

Object::Object() 
{
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif

	//for rendering
	_stepSize = 0.001f;
}

Object::~Object()
{
}

void Object::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
}

void Object::Resize(GLFWwindow* windowHandle, int x, int y)
{
	_arcball.SetWidthHeight(x, y);
}

void Object::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{
	double mouseX, mouseY;

	if(button == GLFW_MOUSE_BUTTON_LEFT) {
		::glfwGetCursorPos(window, &mouseX, &mouseY);
		if(action == GLFW_PRESS) {
			_isLeftKeyPressed = true;
			_arcball.StartRotation(mouseX, mouseY);
		}
		else if(action == GLFW_RELEASE) {
			_isLeftKeyPressed = false;
			_arcball.StopRotation();
		}
	}
	else if(button == GLFW_MOUSE_BUTTON_MIDDLE) {
		_isMiddleKeyPressed = (action == GLFW_PRESS);
	}
	else if(button == GLFW_MOUSE_BUTTON_RIGHT) {
#if 0
		if (action == GLFW_PRESS) {
			_isRightKeyPressed = true;
			::glfwGetCursorPos(window, &mouseX, &mouseY);
			_arcball.StartZooming(mouseX, mouseY);
		}
		else if(action ==GLFW_RELEASE) {
			_isRightKeyPressed = false;
			_arcball.StopZooming();
		}
#endif
	}
}

void Object::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	if(_isLeftKeyPressed && _isCtrlPressed) {
	}
	else if(_isLeftKeyPressed) {
		_arcball.UpdateRotation(nx, ny);
	}
#if 0
	else if(_isRightKeyPressed) {
		_arcball.UpdateZooming(nx, ny);
	}
#endif
}

void Object::MouseScroll(GLFWwindow *window, double nx, double ny)
{
	_arcball.StartZooming(0, 0);
	_arcball.UpdateZooming(-ny, nx);
	_arcball.StopZooming();
}

void Object::Reset() {
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	int width, height;
	glfwGetWindowSize(_windowHandle, &width, &height);
	_winX = width;
	_winY = height;
	_arcball.SetWidthHeight(width, height);

#if 1
	_depth = 3.8;
	_rotX = 10;
	_rotY = -20;


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0, 0, -_depth);
	glRotatef(_rotX, 1, 0, 0);
	glRotatef(_rotY, 0, 1, 0);
#endif
}

void Object::Show() 
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	GLfloat vertices[][3] = { { -1.0f, -1.0f, -1.0f },
		{ 1.0f, -1.0f, -1.0f}, { 1.0f, 1.0f, -1.0f },
		{ -1.0f, 1.0f, -1.0f }, { -1.0f, -1.0f, 1.0f },
		{ 1.0f, -1.0f, 1.0f }, { 1.0f, 1.0f, 1.0f },
		{ -1.0f, 1.0f, 1.0f } };              // Vertices
	int faces[][4] = { { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 4, 5, 6, 7 },
		{ 0, 4, 7, 3 }, { 0, 1, 5, 4 }, { 0, 3, 2, 1 } };

	GLfloat colors[][3] = { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f },
		{ 1.0f, 1.0f, 0.0f }, { 0.0f, 0.5f, 0.5f },
		{ 0.5f, 0.0f, 0.5f }, { 0.5f, 0.5f, 0.0f } };

	for (int i = 0; i < 6; i++) {		
		Eigen::Vector3f nm;
		Eigen::Vector3f v0(vertices[faces[i][0]][0], vertices[faces[i][0]][1], vertices[faces[i][0]][2]);
		Eigen::Vector3f v1(vertices[faces[i][1]][0], vertices[faces[i][1]][1], vertices[faces[i][1]][2]);
		Eigen::Vector3f v2(vertices[faces[i][2]][0], vertices[faces[i][2]][1], vertices[faces[i][2]][2]);
		Eigen::Vector3f l1 = v1-v0;
		Eigen::Vector3f l2 = v2-v0;
		nm = l1.cross(l2);
		nm.normalize();

		glNormal3f(nm(0),nm(1),nm(2));
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < 4; j++)
			glVertex3fv(vertices[faces[i][j]]);
		glEnd();
	}
}
