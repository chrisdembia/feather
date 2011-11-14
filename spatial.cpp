#include <armadillo>
#include <cmath>
#include "spatial.h"

using namespace arma;

mat rx( double theta)
{
	double c = cos(theta);
	double s = sin(theta);
	
	mat RX(3,3);
	
	RX << 1.0 << 0.0 << 0.0 << endr
	   << 0.0 <<  c  <<  s  << endr
	   << 0.0 << -s  <<  c  << endr;
	
	return RX;
}

mat ry( double theta)
{
	double c = cos(theta);
	double s = sin(theta);
	
	mat RY(3,3);
	
	RY <<  c  << 0.0 << -s  << endr
	   << 0.0 << 1.0 << 0.0 << endr
	   <<  s  << 0.0 <<  c  << endr;
	
	return RY; 
}

mat rz( double theta)
{
	double c = cos(theta);
	double s = sin(theta);
	
	mat RY(3,3);
	
	RY <<  c  <<  s  << 0.0 << endr
	   << -s  <<  c  << 0.0 << endr
	   << 0.0 << 0.0 << 1.0 << endr;
	
	return RY; 
}

mat crmp( vec v)
{
	if ( v.n_elem != 3 )
		cout << "function crmp: v.n_elem != 3 " << v << endl;
		
	mat V(3,3);
	
	V << 0.0   << -v(2) <<  v(1) << endr
	  <<  v(2) <<  0.0  << -v(0) << endr
	  << -v(1) <<  v(0) <<  0.0  << endr;
	  
	return V;
}


mat rotx( double theta)
{
	mat E = rx(theta);
	
	mat X = zeros<mat>(6,6);
	
	X( span(0,2), span(0,2) ) = E;
	X( span(3,5), span(3,5) ) = E;
	
	return X;
}

mat roty( double theta)
{
	mat E = ry(theta);
	
	mat X = zeros<mat>(6,6);
	
	X( span(0,2), span(0,2) ) = E;
	X( span(3,5), span(3,5) ) = E;
	
	return X;
}

mat rotz( double theta)
{
	mat E = rz(theta);
	
	mat X = zeros<mat>(6,6);
	
	X( span(0,2), span(0,2) ) = E;
	X( span(3,5), span(3,5) ) = E;
	
	return X;
}

mat xlt( vec r)
{
	if ( r.n_elem != 3 )
		cout << "function xlt: r.n_elem != 3 " << r << endl;
		
	mat X = eye<mat>(6,6);
	
	X( span(3,5), span(0,2) ) = -crmp(r);
	
	return X;
}

mat crm( vec v)
{
	if ( v.n_elem != 6 )
		cout << "function crm: v.n_elem != 6 " << v << endl;
		
	mat V = zeros<mat>(6,6);
	
	V( span(0,2), span(0,2) ) = crmp( v.rows(0,2) );
	V( span(3,5), span(3,5) ) = crmp( v.rows(0,2) );
	V( span(3,5), span(0,2) ) = crmp( v.rows(3,5) );
	
	return V;
}

mat crf( vec v)
{
	if ( v.n_elem != 6 )
		cout << "function crf: v.n_elem != 6 " << v << endl;
		
	return -trans(crm(v));
}

mat mcI( double m, vec c, mat Ic)
{
	if ( c.n_elem != 3 )
		cout << "function mcI: c.n_elem != 3 " << c << endl;
	if ( Ic.n_elem != 9 )
		cout << "function mcI: Ic.n_elem != 9 " << Ic << endl;
		
	mat I(6,6);
	
	I( span(0,2), span(0,2) ) = Ic - m * crmp(c) * crmp(c);
	I( span(0,2), span(3,5) ) =  m * crmp(c);
	I( span(3,5), span(0,2) ) = -m * crmp(c);
	I( span(3,5), span(3,5) ) =  m * eye<mat>(3,3);
	
	return I;
}

vec XtoV( mat X)
{
	if ( X.n_elem != 36 )
		cout << "function XtoV X.n_elem != 36 " << X << endl;
		
	vec v(6);
	v(0) = X(1,2) - X(2,1);
	v(1) = X(2,0) - X(0,2);
	v(2) = X(0,1) - X(1,0);
	v(3) = X(4,2) - X(5,1);
	v(4) = X(5,0) - X(3,2);
	v(5) = X(3,1) - X(4,0);
	
	v = 0.5*v;
	
	return v;
}

