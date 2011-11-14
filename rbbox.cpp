#include <GL/gl.h>
#include <GL/glu.h>
#include <armadillo>
#include "spatial.h"

#include "rb.h"
#include "rbbox.h"

using namespace std;
using namespace arma;

void RBbox::setW( double w)
{
	W = w;
}

void RBbox::setH( double h)
{
	H = h;
}

void RBbox::setD( double d)
{
	D = d;
}

double RBbox::getW()
{
	return W;
}

double RBbox::getH()
{
	return H;
}

double RBbox::getD()
{
	return D;
}

void RBbox::calcCOM()
{
	vec com(3);
	com(0) = W/2;
	com(1) = H/2;
	com(2) = D/2;
	
	RB::setCOM( com);
}

void RBbox::calcMCI()
{
	vec I(3);
	I(0) = 1/12*getMass()*( pow(H,2) + pow(D,2) );
	I(1) = 1/12*getMass()*( pow(W,2) + pow(D,2) );
	I(2) = 1/12*getMass()*( pow(W,2) + pow(H,2) );

	RB::setMCI( mcI( getMass(), getCOM(), diagmat(I)) );
}

void RBbox::draw( vec pos, vec axis)
{ // color?
//	glLoadIdentity();
	glTranslatef( pos(0), pos(1), pos(2));
	glRotatef( 180/math::pi() * sqrt( dot( axis, axis) ) , axis(0), axis(1), axis(2) );

	double w = getW();
	double h = getH()/2;
	double d = getD()/2;
	
	glColor3f( 1.0f, 0.0f, 0.0f);
	glBegin(GL_QUADS);
		// face 1
		glVertex3f( 0.0f, -h  , -d  );
		glVertex3f(  w  , -h  , -d  );
		glVertex3f(  w  ,  h  , -d  );
		glVertex3f( 0.0f,  h  , -d  );		
	glEnd();
	
	glBegin(GL_QUADS);	
		// face 2
		glVertex3f( 0.0f, -h  ,  d  );
		glVertex3f(  w  , -h  ,  d  );
		glVertex3f(  w  ,  h  ,  d  );
		glVertex3f( 0.0f,  h  ,  d  );		
	glEnd();
	
	glColor3f( 0.0f, 0.8f, 0.0f);
	glBegin(GL_QUADS);
		// face 3
		glVertex3f( 0.0f, -h  , -d  );
		glVertex3f(  w  , -h  , -d  );
		glVertex3f(  w  , -h  ,  d  );
		glVertex3f( 0.0f, -h  ,  d  );

	glEnd();
	
	glBegin(GL_QUADS);
		// face 4
		glVertex3f( 0.0f,  h  , -d  );
		glVertex3f(  w  ,  h  , -d  );
		glVertex3f(  w  ,  h  ,  d  );
		glVertex3f( 0.0f,  h  ,  d  );
	glEnd();
	
	glColor3f( 0.0f, 0.0f, 1.0f);
	glBegin(GL_QUADS);
		// face 5
		glVertex3f( 0.0f, -h  , -d  );
		glVertex3f( 0.0f,  h  , -d  );
		glVertex3f( 0.0f,  h  ,  d  );
		glVertex3f( 0.0f, -h  ,  d  );
	glEnd();
	
	glBegin(GL_QUADS);
		// face 6
		glVertex3f(  w  , -h  , -d  );
		glVertex3f(  w  ,  h  , -d  );
		glVertex3f(  w  ,  h  ,  d  );		
		glVertex3f(  w  , -h  ,  d  );
	glEnd();
}


