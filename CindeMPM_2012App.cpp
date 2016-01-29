#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Vbo.h"
#include "cinder/Camera.h"
#include "cinder/gl/GLee.h"
#include "cinder/params/Params.h"
#include "cinder/ImageIo.h"
#include "cinder/Utilities.h"

#include <sstream>

#include <vector>
#include "FF_simulator.h"


using namespace ci;
using namespace ci::app;
using namespace std;

class CindeMPM_2012App : public AppNative {
	simulator s;
	int n;
	GLfloat* vertices;
	CameraPersp m_camera;
	Vec3f mEye, mCenter,mUp;

  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
	void prepareSettings(Settings *settings);
};

void CindeMPM_2012App::setup()
{
	//s.initializeGrid(400,200);
	//s.addParticles();
	//n = s.particles.size();
	//s=new simulator();
	s.initialize_explosion();
	gl::VboMesh::Layout layout;
	layout.setDynamicPositions();
	
	vertices = new GLfloat[s.p_particle_data->particle_list.size()*3];
	m_camera.setPerspective(45.0f, 1.5f, 0.01,100);
	mEye=cinder::Vec3f(-2, 2.8, 4.5)+cinder::Vec3f(-2.75, 2.05, 3.75).normalized()*-2;
	mCenter=cinder::Vec3f(0.75, 0.75, 0.75);
	mUp=cinder::Vec3f(0.27,0.93,-0.21);


}

void CindeMPM_2012App::mouseDown( MouseEvent event )
{
}

void CindeMPM_2012App::update()
{
	//cout<<getHomeDirectory()<<endl;
	writeImage( getHomeDirectory() / "cinder" / "saveImage_" / ( toString( s.p_particle_data->current_step ) + ".png" ), copyWindowSurface() ); 
	int mm=s.p_particle_data->num_particle;
	s.update();
	//mEye = Vec3f( 0.0f, 0.0f, mCameraDistance );
	m_camera.lookAt( mEye, mCenter, mUp );
	gl::setMatrices( m_camera );
	//gl::rotate( mSceneRotation );
}

void CindeMPM_2012App::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) );
	//gl::enableDepthRead();
	//gl::enableDepthWrite();
	gl::color(1,.5,1);
	vector<particle*>& particles = s.p_particle_data->particle_list;
	int nParticles = particles.size();
	float* vi = vertices;
	cinder::Vec3f minp(1000, 1000, 1000);
	cinder::Vec3f maxp(-1000, -1000, -1000);
	for (int i = 0; i < nParticles; i++) {
		particle* p = particles[i];
		*(vi++) = p->position_t.x;
		*(vi++) = p->position_t.y;
		*(vi++) = p->position_t.z;

	}
	cout<<"step:"<<s.p_particle_data->current_step<<endl;
	cout<<"time:"<<s.p_particle_data->current_time<<endl;
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glPointSize(1);

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertices);
	
	//glDrawArrays(GL_LINES, 0, nParticles/2);
	glDrawArrays(GL_POINTS,0,nParticles);

	gl::color(1,1,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	   glVertex3f(0 , 0, 0);

	glEnd();
	glDisableClientState(GL_VERTEX_ARRAY);

	//mParams.draw();
}

void CindeMPM_2012App::prepareSettings( Settings *settings ) {
	settings->setWindowSize( 1200, 800 );
	settings->setFrameRate(60.0f);
	settings->enableConsoleWindow(true);
}

CINDER_APP_NATIVE( CindeMPM_2012App, RendererGl )
