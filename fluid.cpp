#include "fluid.h"

Fluid::Fluid()
{
	diffusion = 0.00001f;
	viscosity = 0.000000f;
	buoyancy  = 4.0f;
	vc_eps    = 5.0f;
	_dt = DT;
	_isLightSelected = false;
	_isRendering = true;
	_isDrawSliceOutline = false;

	for (int i=0; i<10; i++)
		ClearBuffer(_buffers[i]);

	int i=0;
	_density=_buffers[i++]; _densityTmp=_buffers[i++];
	_velX=_buffers[i++]; _velXTmp=_buffers[i++];
	_velY=_buffers[i++]; _velYTmp=_buffers[i++];
	_velZ=_buffers[i++]; _velZTmp=_buffers[i++];

	ClearSources();

	_lightPos[0] = -1.2f;
	_lightPos[1] = 0.2f;
	_lightPos[2] = 1.2f;

	_textureData = NULL;
	_isDrawSliceOutline = false;
	_isRendering = true;

	_RES = RES;
	_N = _RES - 2;
	_SIZE = _RES*_RES*_RES;
	_volumeData = _density;

	// cube vertices
	GLfloat cv[][3] = {
		{1.0f, 1.0f, 1.0f}, {-1.0f, 1.0f, 1.0f}, {-1.0f, -1.0f, 1.0f}, {1.0f, -1.0f, 1.0f},
		{1.0f, 1.0f, -1.0f}, {-1.0f, 1.0f, -1.0f}, {-1.0f, -1.0f, -1.0f}, {1.0f, -1.0f, -1.0f}
	};

	// cube edges have the form edges[n][0][xyz] + t*edges[n][1][xyz]
	GLfloat ce[12][2][3] = {
		{{1.0f, 1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f, 1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},

		{{1.0f, -1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},
		{{1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},

		{{-1.0f, 1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, 1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}}
	};

	memcpy(_cubeVertices, cv, sizeof(_cubeVertices) );
	memcpy(_cubeEdges, ce, sizeof(_cubeEdges) );

	InitGL();
	
	_lightDir = -_lightPos;
	GenerateRayTemplate(_RES);
}

class Convexcomp
{
	private:
		const Eigen::Vector3f &p0, &up;
	public:
		Convexcomp(const Eigen::Vector3f& p0, const Eigen::Vector3f& up) : p0(p0), up(up) {}

		bool operator()(const Eigen::Vector3f& a, const Eigen::Vector3f& b) const
		{
			Eigen::Vector3f va = a-p0, vb = b-p0;
			//return dot(up, cross(va, vb)) >= 0;
			return up.dot(va.cross(vb)) >= 0;
		}
};

void Fluid::InitGL()
{
	glEnable(GL_TEXTURE_3D);
	glDisable(GL_DEPTH_TEST);
	glCullFace(GL_FRONT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	glGenTextures(2, &_hTexture);

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glBindTexture(GL_TEXTURE_3D, _hTexture);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R,GL_CLAMP);

}

const float* Fluid::GetDensity()
{
	return _density;
}

bool Fluid::LightSelected(double mouseX, double mouseY)
{
	GLdouble mvMatrix[16], projMatrix[16];
	GLint viewportMatrix[4];
	double objX = 0, objY = 0, objZ = 0;
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewportMatrix);

	//point on 'near' plane
	gluUnProject(mouseX, _winY-mouseY, 0.0,  mvMatrix,  projMatrix,
			viewportMatrix, &objX, &objY, &objZ);
	Eigen::Vector3f ptNear(objX,objY, objZ); 

	//point on 'far' plane
	gluUnProject(mouseX, _winY-mouseY, 1.0,  mvMatrix,  projMatrix,
			viewportMatrix, &objX, &objY, &objZ);
	Eigen::Vector3f ptFar(objX,objY, objZ); 

#if 1
	//calculate distance between point and line
	Eigen::Vector3f crossProduct = (_lightPos-ptNear).cross(_lightPos-ptFar);
	float dist = crossProduct.norm() /
		(ptFar-ptNear).norm();
#else
	//another method: calculate distance between point and line
	Eigen::Vector3f lineDir = (ptFar-ptNear);
	lineDir.normalize();
	Eigen::Vector3f pointDir = (_lightPos-ptNear);
	float t = pointDir.dot(lineDir);
	Eigen::Vector3f projection = ptNear + (lineDir * t);
	float dist = (projection-_lightPos).norm();
#endif
	if (dist < 0.1)
		return true;

	return false;
}

void Fluid::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{
	double mouseX, mouseY;
	if(button == GLFW_MOUSE_BUTTON_LEFT) {
		::glfwGetCursorPos(window, &mouseX, &mouseY);
		if(action == GLFW_PRESS) {
			_isLeftKeyPressed = true;

			if (LightSelected(mouseX, mouseY)) {
				std::cout << "light selected" << std::endl;
				_isLightSelected = true;
				_arcball.StartDragging(mouseX, mouseY);
			}
			_arcball.StartRotation(mouseX, mouseY);
		}
		else if(action == GLFW_RELEASE) {
			_isLeftKeyPressed = false;
			_isLightSelected = false;
			_arcball.StopRotation();
			_arcball.StopDragging();
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

void Fluid::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	if(_isLeftKeyPressed && _isCtrlPressed) {
	}
	else if(_isLeftKeyPressed && _isLightSelected) {
		_lightPos += _arcball.UpdateDragging(nx, ny);
		_lightDir = -_lightPos;
		GenerateRayTemplate(_RES);
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

void Fluid::Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch(key) {
			case GLFW_KEY_S:		// Rendering on/off
				_isRendering = !_isRendering;
				break;
			case GLFW_KEY_W:		//wireframe on/off
				_isDrawSliceOutline = !_isDrawSliceOutline;
				break;
		}
	}
}

Fluid::~Fluid()
{
		if (_textureData)
		free(_textureData);
}

#define BOUNDARY
#undef BOUNDARY
void Fluid::EnforceBoundary(int b, float* quantity)
{
#ifdef BOUNDARY
	int i, j;
	for (i=1; i<=N; i++)
	{
		for (j=1; j<=N; j++) {
			quantity[_I(0,i,j)]    = (b==1) ? -quantity[_I(1,i,j)] : quantity[_I(1,i,j)];
			quantity[_I(N+1,i,j)]  = (b==1) ? -quantity[_I(N,i,j)] : quantity[_I(N,i,j)];
			quantity[_I(i,0,j)]    = (b==2) ? -quantity[_I(i,1,j)] : quantity[_I(i,1,j)];
			quantity[_I(i,N+1,j)]  = (b==2) ? -quantity[_I(i,N,j)] : quantity[_I(i,N,j)];
			quantity[_I(i,j,0)]    = (b==3) ? -quantity[_I(i,j,1)] : quantity[_I(i,j,1)];
			quantity[_I(i,j,N+1)]  = (b==3) ? -quantity[_I(i,j,N)] : quantity[_I(i,j,N)];
		}
	}

	quantity[_I(0,0,0)]       = (quantity[_I(1,0,0)]    + quantity[_I(0,1,0)]     + quantity[_I(0,0,1)])    / 3;
	quantity[_I(0,N+1,0)]     = (quantity[_I(1,N+1,0)]  + quantity[_I(0,N,0)]     + quantity[_I(0,N+1,1)])  / 3;
	quantity[_I(N+1,0,0)]     = (quantity[_I(N,0,0)]    + quantity[_I(N+1,1,0)]   + quantity[_I(N+1,0,1)])  / 3;
	quantity[_I(N+1,N+1,0)]   = (quantity[_I(N,N+1,0)]  + quantity[_I(N+1,N,0)]   + quantity[_I(N+1,N+1,1)])/ 3;
	quantity[_I(0,0,N+1)]     = (quantity[_I(1,0,N+1)]  + quantity[_I(0,1,N+1)]   + quantity[_I(0,0,N)])    / 3;
	quantity[_I(0,N+1,N+1)]   = (quantity[_I(1,N+1,N+1)]+ quantity[_I(0,N,N+1)]   + quantity[_I(0,N+1,N)])  / 3;
	quantity[_I(N+1,0,N+1)]   = (quantity[_I(N,0,N+1)]  + quantity[_I(N+1,1,N+1)] + quantity[_I(N+1,0,N)])  / 3;
	quantity[_I(N+1,N+1,N+1)] = (quantity[_I(N,N+1,N+1)]+ quantity[_I(N+1,N,N+1)] + quantity[_I(N+1,N+1,N)])/ 3;
#endif
}

void Fluid::AddSource(float* src, float *dst)
{
	int i;

	for (i=0; i<SIZE; i++)
		dst[i] += src[i]*_dt;
}

void Fluid::AddBuoyancy()
{
	int i;

	for (i=0; i<SIZE; i++)
		_velY[i] += _density[i]*buoyancy*_dt;//FIXME
}

inline void Fluid::Diffuse(int b, float* velTmp, float* vel, float diff)
{
	int l;
	float a=_dt*diff*N*N*N;
	for (l=0; l<20; l++) 
	{
		FOR_ALL_CELL {
			vel[_I(i,j,k)] = (velTmp[_I(i,j,k)] + a*(
						vel[_I(i-1,j,k)]+vel[_I(i+1,j,k)]+
						vel[_I(i,j-1,k)]+vel[_I(i,j+1,k)]+
						vel[_I(i,j,k-1)]+vel[_I(i,j,k+1)]))/(1+6*a);
		} END_FOR
		EnforceBoundary(b,vel);
	}
}

inline void Fluid::Advect(int b, float* quantityTmp, float* quantity, float* velX, float* velY, float* velZ)
{
	int i0, j0, k0, i1, j1, k1;
	float sx0, sx1, sy0, sy1, sz0, sz1, v0, v1;
	float xx, yy, zz, dt0;
	dt0 = _dt*N;
	FOR_ALL_CELL {
		xx = i - dt0*velX[_I(i,j,k)];
		yy = j - dt0*velY[_I(i,j,k)];
		zz = k - dt0*velZ[_I(i,j,k)];
		if (xx<0.5) xx=0.5f; if (xx>N+0.5) xx=N+0.5f; i0=(int)xx; i1=i0+1;
		if (yy<0.5) yy=0.5f; if (yy>N+0.5) yy=N+0.5f; j0=(int)yy; j1=j0+1;
		if (zz<0.5) zz=0.5f; if (zz>N+0.5) zz=N+0.5f; k0=(int)zz; k1=k0+1;
		sx1 = xx-i0; sx0 = 1-sx1;
		sy1 = yy-j0; sy0 = 1-sy1;
		sz1 = zz-k0; sz0 = 1-sz1;
		v0 = sx0 * ( sy0*quantityTmp[_I(i0,j0,k0)] + sy1*quantityTmp[_I(i0,j1,k0)] ) +
			sx1 * ( sy0*quantityTmp[_I(i1,j0,k0)] + sy1*quantityTmp[_I(i1,j1,k0)] );
		v1 = sx0 * ( sy0*quantityTmp[_I(i0,j0,k1)] + sy1*quantityTmp[_I(i0,j1,k1)] ) +
			sx1 * ( sy0*quantityTmp[_I(i1,j0,k1)] + sy1*quantityTmp[_I(i1,j1,k1)] );
		quantity[_I(i,j,k)] = sz0*v0 + sz1*v1;
	} END_FOR
	EnforceBoundary(b,_density);
}

void Fluid::Project(void)
{
	float* p = _velXTmp;	
	float* div = _velYTmp;	// temporary buffers, use old velocity buffers
	int l;
	float h;
	h = 1.0f/N;
	FOR_ALL_CELL {
		div[_I(i,j,k)] = -h*(
				_velX[_I(i+1,j,k)]-_velX[_I(i-1,j,k)]+
				_velY[_I(i,j+1,k)]-_velY[_I(i,j-1,k)]+
				_velZ[_I(i,j,k+1)]-_velZ[_I(i,j,k-1)])/3;
		p[_I(i,j,k)] = 0;
	} END_FOR
	EnforceBoundary(0,div); EnforceBoundary(0,p);
	for (l=0; l<20; l++) 
	{
		FOR_ALL_CELL {
			p[_I(i,j,k)] = (div[_I(i,j,k)]+
					p[_I(i-1,j,k)]+p[_I(i+1,j,k)]+
					p[_I(i,j-1,k)]+p[_I(i,j+1,k)]+
					p[_I(i,j,k-1)]+p[_I(i,j,k+1)])/6;
		} END_FOR
		EnforceBoundary(0,p);
	}
	FOR_ALL_CELL {
		_velX[_I(i,j,k)] -= (p[_I(i+1,j,k)]-p[_I(i-1,j,k)])/3/h;
		_velY[_I(i,j,k)] -= (p[_I(i,j+1,k)]-p[_I(i,j-1,k)])/3/h;
		_velZ[_I(i,j,k)] -= (p[_I(i,j,k+1)]-p[_I(i,j,k-1)])/3/h;
	} END_FOR
	EnforceBoundary(1,_velX); EnforceBoundary(2,_velY);
}

void Fluid::VorticityConfinement()
{
	int ijk;
	//temp buffers
	float *curlX = _velXTmp, *curlY = _velYTmp, *curlZ=_velZTmp, *curl=_densityTmp;
	float dt0 = _dt * vc_eps;

	FOR_ALL_CELL {
		ijk = _I(i,j,k);
		// curlx = dw/dy - dv/dz
		curlX[ijk] = (_velZ[_I(i,j+1,k)] - _velZ[_I(i,j-1,k)]) * 0.5f -
			(_velY[_I(i,j,k+1)] - _velY[_I(i,j,k-1)]) * 0.5f;

		// curly = du/dz - dw/dx
		curlY[ijk] = (_velX[_I(i,j,k+1)] - _velX[_I(i,j,k-1)]) * 0.5f -
			(_velZ[_I(i+1,j,k)] - _velZ[_I(i-1,j,k)]) * 0.5f;

		// curlz = dv/dx - du/dy
		curlZ[ijk] = (_velY[_I(i+1,j,k)] - _velY[_I(i-1,j,k)]) * 0.5f -
			(_velX[_I(i,j+1,k)] - _velX[_I(i,j-1,k)]) * 0.5f;

		// curl = |curl|
		curl[ijk] = sqrtf(curlX[ijk]*curlX[ijk] +
				curlY[ijk]*curlY[ijk] +
				curlZ[ijk]*curlZ[ijk]);
	} END_FOR

	FOR_ALL_CELL {
		ijk = _I(i,j,k);
		float nX = (curl[_I(i+1,j,k)] - curl[_I(i-1,j,k)]) * 0.5f;
		float nY = (curl[_I(i,j+1,k)] - curl[_I(i,j-1,k)]) * 0.5f;
		float nZ = (curl[_I(i,j,k+1)] - curl[_I(i,j,k-1)]) * 0.5f;
		float len1 = 1.0f/(sqrtf(nX*nX+nY*nY+nZ*nZ)+0.0000001f);
		nX *= len1;
		nY *= len1;
		nZ *= len1;
		_velX[ijk] += (nY*curlZ[ijk] - nZ*curlY[ijk]) * dt0;
		_velY[ijk] += (nZ*curlX[ijk] - nX*curlZ[ijk]) * dt0;
		_velZ[ijk] += (nX*curlY[ijk] - nY*curlX[ijk]) * dt0;
	} END_FOR
}

#define DIFFUSE
#define ADVECT

void Fluid::VelocityStep()
{
#if 1
	AddSource(su, _velX);
	AddSource(sv, _velY);
	AddSource(sw, _velZ);
#endif
	AddBuoyancy();
	VorticityConfinement();

#ifdef DIFFUSE
	std::swap(_velXTmp, _velX); std::swap(_velYTmp, _velY); std::swap(_velZTmp, _velZ);
	Diffuse(1, _velXTmp, _velX, viscosity);
	Diffuse(2, _velYTmp, _velY, viscosity);
	Diffuse(3, _velZTmp, _velZ, viscosity);
	Project();
#endif
#ifdef ADVECT
	std::swap(_velXTmp, _velX); std::swap(_velYTmp, _velY); std::swap(_velZTmp, _velZ);
	Advect(1, _velXTmp, _velX, _velXTmp, _velYTmp, _velZTmp);
	Advect(2, _velYTmp, _velY, _velXTmp, _velYTmp, _velZTmp);
	Advect(3, _velZTmp, _velZ, _velXTmp, _velYTmp, _velZTmp);
	Project();
#endif
}

void Fluid::DensityStep()
{
#if 1
	AddSource(sd, _density);
#endif
#ifdef DIFFUSE
	std::swap(_densityTmp, _density);
	Diffuse(0, _densityTmp, _density, diffusion);
#endif
#ifdef ADVECT
	std::swap(_densityTmp, _density);
	Advect(0, _densityTmp, _density, _velX, _velY, _velZ);

#if 1
	//decrease density
	FOR_ALL_CELL {
		_density[_I(k,j,i)] -= 0.001;
		if(_density[_I(k,j,i)] < 0)
			_density[_I(k,j,i)] = 0;
	} END_FOR
#endif

#endif
}

void Fluid::GenerateSmoke()
{
	const int centerY = RES/4;
	const int centerZ = RES/2;
	float dens = (rand()%1000)/1000.0f;
	for (int i = 1 ; i <= N ; i++) {
		for (int j = 1 ; j <= N ; j++) {
			if (hypot(i-centerY, j-centerZ) < RES/10) {
				this->_density[_I(5,i,j)] = dens;
				this->_velX[_I(5,i,j)] = 2.0f;
			}
		}
	}

}

void Fluid::SimulateStep()
{
	GenerateSmoke();

	VelocityStep();
	DensityStep();

}

void Fluid::Render(void)
{
	GLdouble mvMatrix[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

	glDisable(GL_DEPTH_TEST);	//Important in this rendering

	DrawCube();
	DrawLight();
	if(_isRendering) {
		//core rendering parts
		DrawSlices(mvMatrix);
	}
	else 
		DrawVolumeData();

	glEnable(GL_DEPTH_TEST);
}


void Fluid::DrawLight()
{
	
	glPushAttrib(GL_POINT_BIT | GL_CURRENT_BIT);	//save current state
	//draw light
	glPointSize(13.0f);
	glBegin(GL_POINTS);
	glColor4f(0.0f, 1.0f, 1.0f, 1.0f);
	glVertex3f(_lightPos[0], _lightPos[1], _lightPos[2]);
	glEnd();
	glPopAttrib();		//restore painting state
}

void Fluid::DrawVolumeData()
{
	glBegin(GL_POINTS);
	FOR_ALL_CELL {
		if(!ALMOST_EQUAL(_volumeData[_I(i,j,k)], 0) ) {
			glVertex3f(((float)i/_N-0.5)*2, ((float)j/_N-0.5)*2, ((float)k/_N-0.5)*2 );
		}
	} END_FOR
	glEnd();

	glPointSize(13.0f);
	glBegin(GL_POINTS);
	glColor4f(0.0f, 1.0f, 1.0f, 1.0f);
	glVertex3f(_lightPos[0], _lightPos[1], _lightPos[2]);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glEnd();
	glPointSize(1.0f);

}

void Fluid::DrawCube(void)
{
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_FRAGMENT_PROGRAM_ARB);

	glEnable(GL_CULL_FACE);
	float coeff = 0.5;
	glColor4f(0.9f*coeff, 0.8f*coeff, 0.4f*coeff, 1.0f);
	GLfloat (*cv)[3] = _cubeVertices;

	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]); glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[4]); glVertex3fv(cv[5]); glVertex3fv(cv[1]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[3]); glVertex3fv(cv[7]); glVertex3fv(cv[4]);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3fv(cv[7]); glVertex3fv(cv[6]); glVertex3fv(cv[5]); glVertex3fv(cv[4]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]); glVertex3fv(cv[7]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]); glVertex3fv(cv[6]); glVertex3fv(cv[2]);
	glEnd();

	glDisable(GL_CULL_FACE);

	glBegin(GL_LINES);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]);
	glVertex3fv(cv[1]); glVertex3fv(cv[2]);
	glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glVertex3fv(cv[3]); glVertex3fv(cv[0]);

	glVertex3fv(cv[4]); glVertex3fv(cv[5]);
	glVertex3fv(cv[5]); glVertex3fv(cv[6]);
	glVertex3fv(cv[6]); glVertex3fv(cv[7]);
	glVertex3fv(cv[7]); glVertex3fv(cv[4]);

	glVertex3fv(cv[0]); glVertex3fv(cv[4]);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]);
	glVertex3fv(cv[3]); glVertex3fv(cv[7]);
	glEnd();

}

void Fluid::DrawSlices(GLdouble mvMatrix[16])
{
	int i;
	Eigen::Vector3f viewdir(-mvMatrix[2], -mvMatrix[6], -mvMatrix[10]);	//FIXME

	viewdir.normalize();
	// find cube vertex that is closest to the viewer
	GLfloat (*cv)[3] = _cubeVertices;
	for (i=0; i<8; i++) {
		float x = cv[i][0] + viewdir[0];
		float y = cv[i][1] + viewdir[1];
		float z = cv[i][2] + viewdir[2];
		if ((x>=-1.0f)&&(x<=1.0f)
				&&(y>=-1.0f)&&(y<=1.0f)
				&&(z>=-1.0f)&&(z<=1.0f))
		{
			break;
		}
	}
	assert(i != 8);

	// our slices are defined by the plane equation A*x + B*y +C*z + D = 0
	// (a,b,c), the plane normal, are given by viewdir
	// d is the parameter along the view direction. the first d is given by
	// inserting previously found vertex into the plane equation
	float d0 = -(viewdir[0]*cv[i][0] + viewdir[1]*cv[i][1] + viewdir[2]*cv[i][2]);
	float dStep = 2*d0/SLICE_NUM;
	int n = 0;
	for (float d = -d0; d < d0; d += dStep) {
		// IntersectEdges returns the intersection points of all cube edges with
		// the given plane that lie within the cube
		std::vector<Eigen::Vector3f> pt = IntersectEdges(viewdir[0], viewdir[1], viewdir[2], d);

		if (pt.size() > 2) {
			// sort points to get a convex polygon
			std::sort(pt.begin()+1, pt.end(), Convexcomp(pt[0], viewdir));

			glEnable(GL_TEXTURE_3D);
			glBegin(GL_POLYGON);
			for (i=0; i<pt.size(); i++){
				glColor3f(1.0, 1.0, 1.0);
				glTexCoord3d((pt[i][0]+1.0)/2.0, (pt[i][1]+1)/2.0, (pt[i][2]+1.0)/2.0);//FIXME
				glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
			}
			glEnd();

			if (_isDrawSliceOutline)
			{
				glDisable(GL_TEXTURE_3D);
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glBegin(GL_POLYGON);
				for (i=0; i<pt.size(); i++) {
					glColor3f(0.0, 0.0, 1.0);
					glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
				}
				glEnd();
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
		}
		n++;
	}

	glDisable(GL_TEXTURE_3D);
}

std::vector<Eigen::Vector3f> Fluid::IntersectEdges(float A, float B, float C, float D)
{
	float t;
	Eigen::Vector3f p;
	std::vector<Eigen::Vector3f> res;
	GLfloat (*edges)[2][3] = _cubeEdges;
	

	for (int i=0; i<12; i++) {
		t = -(A*edges[i][0][0] + B*edges[i][0][1] + C*edges[i][0][2] + D)
			/ (A*edges[i][1][0] + B*edges[i][1][1] + C*edges[i][1][2]);
		if ((t>0)&&(t<2)) {
			p[0] = edges[i][0][0] + edges[i][1][0]*t;
			p[1] = edges[i][0][1] + edges[i][1][1]*t;
			p[2] = edges[i][0][2] + edges[i][1][2]*t;
			res.push_back(p);
		}
	}

	return res;
}


void Fluid::FillTexture()
{
	if (_textureData == NULL) {
		_textureData = new unsigned char[_SIZE*4];
	}
	memset(_textureData, 0, _SIZE*4);

	const float *density = _volumeData;
	unsigned char* intensity = new unsigned char[_SIZE];
	memset(intensity, 0, _SIZE);
	CastLight(_RES, density, intensity);

#if 1
	//FIXME: It is important to beware that texture coordinate
	//is in reverse order  of the simulation coordinate
	for (int i = 0 ; i < _RES ; i++) {
		for (int j = 0 ; j < _RES ; j++) {
			for (int k = 0 ; k < _RES ; k++) {
				//unsigned char c = 200;
				int texIndex = i*_RES*_RES + j*_RES + k;	/*reverse order*/
				int densIndex = k*_RES*_RES + j*_RES + i;
				unsigned char c = intensity[densIndex];
				_textureData[texIndex*4] = c;
				_textureData[texIndex*4+1] = c;
				_textureData[texIndex*4+2] = c;
				_textureData[texIndex*4+3] = (density[densIndex]>0.1f) ? 
					255 : (unsigned char) (density[densIndex]*2550.0f);
			}
		}
	}
#else
	for (int i=0; i<_SIZE; i++) {
		unsigned char c = intensity[i];
		//unsigned char c = 200;
		_textureData[(i<<2)] = c;
		_textureData[(i<<2)+1] = c;
		_textureData[(i<<2)+2] = c;
		_textureData[(i<<2)+3] = (density[i]>0.1f) ? 255 : (unsigned char) (density[i]*2550.0f);
	}
#endif

	delete []intensity;

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, _RES, _RES, _RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, _textureData);
}

void Fluid::GenerateRayTemplate(int edgeLen)
{
	_rayTemplate[0][0] = _rayTemplate[0][2] = _rayTemplate[0][2] = 0;
	float fx = 0.0f, fy = 0.0f, fz = 0.0f;
	int x = 0, y = 0, z = 0;
	float lx = _lightDir[0] + 0.000001f, ly = _lightDir[1] + 0.000001f, lz = _lightDir[2] + 0.000001f;
	int xinc = (lx > 0) ? 1 : -1;
	int yinc = (ly > 0) ? 1 : -1;
	int zinc = (lz > 0) ? 1 : -1;
	float tx, ty, tz;
	int i = 1;
	int len = 0;
	int maxlen = 3*edgeLen*edgeLen;
	while (len <= maxlen)
	{
		// fx + t*lx = (x+1)   ->   t = (x+1-fx)/lx
		tx = (x+xinc-fx)/lx;
		ty = (y+yinc-fy)/ly;
		tz = (z+zinc-fz)/lz;

		if ((tx<=ty)&&(tx<=tz)) {
			_rayTemplate[i][0] = _rayTemplate[i-1][0] + xinc;
			x =+ xinc;
			fx = x;

			if (ALMOST_EQUAL(ty,tx)) {
				_rayTemplate[i][1] = _rayTemplate[i-1][1] + yinc;
				y += yinc;
				fy = y;
			} else {
				_rayTemplate[i][1] = _rayTemplate[i-1][1];
				fy += tx*ly;
			}

			if (ALMOST_EQUAL(tz,tx)) {
				_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_rayTemplate[i][2] = _rayTemplate[i-1][2];
				fz += tx*lz;
			}
		} else if ((ty<tx)&&(ty<=tz)) {
			_rayTemplate[i][0] = _rayTemplate[i-1][0];
			fx += ty*lx;

			_rayTemplate[i][1] = _rayTemplate[i-1][1] + yinc;
			y += yinc;
			fy = y;

			if (ALMOST_EQUAL(tz,ty)) {
				_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_rayTemplate[i][2] = _rayTemplate[i-1][2];
				fz += ty*lz;
			}
		} else {
			assert((tz<tx)&&(tz<ty));
			_rayTemplate[i][0] = _rayTemplate[i-1][0];
			fx += tz*lx;
			_rayTemplate[i][1] = _rayTemplate[i-1][1];
			fy += tz*ly;
			_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
			z += zinc;
			fz = z;
		}

		len = _rayTemplate[i][0]*_rayTemplate[i][0]
			+ _rayTemplate[i][1]*_rayTemplate[i][1]
			+ _rayTemplate[i][2]*_rayTemplate[i][2];
		i++;
	}
}

#define DECAY 0.06f
void Fluid::CastLight(int n /*edgelen*/, const float* dens, unsigned char* intensity)
{

	int i,j;
	int sx = (_lightDir[0]>0) ? 0 : n-1;
	int sy = (_lightDir[1]>0) ? 0 : n-1;
	int sz = (_lightDir[2]>0) ? 0 : n-1;

	float decay = 1.0f/(n*DECAY);

	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			if (!ALMOST_EQUAL(_lightDir[0], 0) )
				LightRay(sx, i, j, n, decay, dens, intensity);
			if (!ALMOST_EQUAL(_lightDir[1], 0) )
				LightRay(i, sy, j, n, decay, dens, intensity);
			if (!ALMOST_EQUAL(_lightDir[2], 0) )
				LightRay(i, j, sz, n, decay, dens, intensity);
		}
}


#define AMBIENT 100
inline void Fluid::LightRay(int x, int y, int z, int n, float decay, const float* dens, unsigned char* intensity)
{
	int xx = x, yy = y, zz = z, i = 0;
	int offset;

	int l = 200;
	float d;

	do {
		offset = ((xx*n) + yy)*n + zz;//FIXME
		if (intensity[offset] > 0)
			intensity[offset] = (unsigned char) ((intensity[offset] + l)*0.5f);
		else
			intensity[offset] = (unsigned char) l;
		d = dens[offset]*255.0f;
		if (l > AMBIENT) {
			l -= d*decay;
			if (l < AMBIENT)
				l = AMBIENT;
		}

		i++;
		xx = x + _rayTemplate[i][0];
		yy = y + _rayTemplate[i][1];
		zz = z + _rayTemplate[i][2];
	} while ((xx>=0)&&(xx<n)&&(yy>=0)&&(yy<n)&&(zz>=0)&&(zz<n));
}

void Fluid::ClearBuffer(float* buf)
{
	for (int i=0; i<SIZE; i++) {
		buf[i] = 0.0f;
	}
}

void Fluid::ClearSources(void)
{
	for (int i=0; i<SIZE; i++) {
		sd[i] = su[i] = sv[i] = sw[i] = 0.0f;
	}
}

void Fluid::Reset() {
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

void Fluid::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
}

void Fluid::Resize(GLFWwindow* windowHandle, int x, int y)
{
	_arcball.SetWidthHeight(x, y);
}

void Fluid::MouseScroll(GLFWwindow *window, double nx, double ny)
{
	_arcball.StartZooming(0, 0);
	_arcball.UpdateZooming(-ny, nx);
	_arcball.StopZooming();
}


void Fluid::Show()
{
	FillTexture();
	Render();
}

