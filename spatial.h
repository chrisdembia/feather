#ifndef SPATIAL_H
#define SPATIAL_H
/*
#include <armadillo>
#include <cmath>
*/
using namespace arma;

mat rx( double theta);

mat ry( double theta);

mat rz( double theta);

mat crmp( vec v);

mat rotx( double theta);

mat roty( double theta);

mat rotz( double theta);

mat xlt( vec r);

mat crm( vec v);

mat crf( vec v);

mat mcI( double m, vec c, mat Ic);

vec XtoV( mat X);

#endif
