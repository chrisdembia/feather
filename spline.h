#ifndef SPLINE_H
#define SPLINE_H

using namespace std;
using namespace arma;

class spline
{
	public:
	spline(){}
	spline( char* fn)
	{
		readViaPts( fn );

		coeffs.set_size(N);

		calcCoeffs();
		
	}
	
	int getNvp();
	int getN();
	vec getT();
	mat getP();
	void readViaPts(char* fname);
	mat calcCoeff( vec t, vec p);
	void calcCoeffs();
	void calcCoeffs( vec t, mat p);
	
	vec eval( double t);
	double eval( int i, double t);
	
	private:
	int Nvp; // number of via points
	int N; // number of splines (N dof of manipulator)
	vec T; // time at via points, length Nvp
	mat P;
	field<mat> coeffs;
	
	void setN(int n);
	void setNvp(int nvp);
	
	int getSplineNum( double &t);
	double getAQ(int n, int i, double t);
};

#endif
