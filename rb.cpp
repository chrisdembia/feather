#include <cstdlib>
#include <GL/gl.h>
#include <GL/glu.h>
#include <armadillo>

#include "rb.h"

using namespace std;
using namespace arma;

void RB::setMass( double m)
{
	mass = m;
}

void RB::setCOM( vec COM)
{
	com = COM;
}

void RB::setMCI( mat mci)
{
	MCI = mci;
}

double RB::getMass()
{
	return mass;
}

vec RB::getCOM()
{
	return com;
}

mat RB::getMCI()
{
	return MCI;
}

