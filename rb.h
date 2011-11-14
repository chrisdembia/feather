#ifndef RB_H
#define RB_H
/*
#include <cstdlib>
#include <GL/gl.h>
#include <GL/glu.h>
#include <armadillo>
#include "SDL.h"
*/

using namespace std;
using namespace arma;

class RB
{
	public:
	RB(){};
	
	void setMass( double m);
	void setCOM( vec COM);
	void setMCI( mat mci);


	double getMass();
	vec getCOM();
	mat getMCI();
	

	virtual double getR(){ return 20; }
	
	void calcCOM();
	virtual void draw( vec x0, vec axis){}

	private:
	double mass;
	vec com;
	mat MCI;

};


#endif
