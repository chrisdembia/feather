#ifndef PROPHET_H
#define PROPHET_H
/*
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <armadillo>
*/
using namespace std;
using namespace arma;
   
int crystalball(double t, const double CQ[], double Cdots[], void *params);
     
mat prophet(model &M, vec tspan, vec ICs, vec Ks, spline &sp1);

void ARMAtoC( vec v, double y[]);
	
vec CtoARMA( double v[], int N);


#endif
