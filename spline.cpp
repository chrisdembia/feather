#include <cstdlib>
#include <armadillo>
#include <algorithm>

#include "spline.h"

using namespace std;
using namespace arma;

int spline::getNvp()
{
	return Nvp;
}

int spline::getN()
{
	return N;
}

vec spline::getT()
{
	return T;
}

mat spline::getP()
{
	return P;
}

vec spline::eval( double t)
{
	vec q(N);

	int n = getSplineNum( t);
	
	for (int i = 0; i < N; i++)
	{
		q(i) = getAQ( n, i, t);
	}
	
	return q;
}

double spline::eval(int i, double t)
{
	int n = getSplineNum( t);
	
	// use proper coefficients
	return getAQ( n, i, t);
}

void spline::setN( int n)
{
	N = n;
}

void spline::setNvp( int nvp)
{
	Nvp = nvp;
}


mat spline::calcCoeff( vec t, vec p)
{
	int np = p.n_elem;
	
	int N = 4 * (np - 1);
	
	mat M = zeros<mat>(N,N);
	
	vec b = zeros<vec>(N);
	
	b(0) = p(0);

	M(0,0) = 1;
	M(0,1) = t(0);
	M(0,2) = t(0)*t(0);
	M(0,3) = t(0)*t(0)*t(0);
	
	M(1,1) = 1;
	M(1,2) = 2*t(0);
	M(1,3) = 3*t(0)*t(0);
	
	M(N-2,N-4) = 1;
	M(N-2,N-3) = t(np-1);
	M(N-2,N-2) = t(np-1)*t(np-1);
	M(N-2,N-1) = t(np-1)*t(np-1)*t(np-1);
	
	M(N-1,N-3) = 1;
	M(N-1,N-2) = 2*t(np-1);
	M(N-1,N-1) = 3*t(np-1)*t(np-1);
	
	b(N-2) = p(np-1);

	// manage center points
	int idx = 1;
	double th;
	for (int i = 1; i < N-4; i+=4 )
	{
		th = t(idx);

		// solution vector
		b(i+1) = p(idx);
		b(i+4) = p(idx);
	
		// position on the left
		M(i+1,i-1) = 1;
		M(i+1,i) = th;
		M(i+1,i+1) = th*th;
		M(i+1,i+2) = th*th*th;
		
		// velocity on left = velocity on right
		M(i+2,i) = 1;
		M(i+2,i+1) = 2*th;
		M(i+2,i+2) = 3*th*th;
		M(i+2,i+4) = -1;
		M(i+2,i+5) = -2*th;
		M(i+2,i+6) = -3*th*th;
		
		// acceleration on the left = acceleration on the right
		M(i+3,i+1) = 2;
		M(i+3,i+2) = 6*th;
		M(i+3,i+5) = -2;
		M(i+3,i+6) = -6*th;
		
		// position on the right
		M(i+4,i+3) = 1;
		M(i+4,i+4) = th;
		M(i+4,i+5) = th*th;
		M(i+4,i+6) = th*th*th;
			
		idx++;
	}

	// solve
	vec a = solve( M, b);

	mat coeff(np-1,4);
	for (int i = 0; i < np-1; i++)
	{
		coeff(i,0) = a(4*i+0);
		coeff(i,1) = a(4*i+1);
		coeff(i,2) = a(4*i+2);
		coeff(i,3) = a(4*i+3);
	}
	
	return coeff;
}

void spline::calcCoeffs()
{

	for (int i = 0; i < N; i++)
	{
		coeffs(i) = calcCoeff( T, P.col(i) );
	}
}

void spline::calcCoeffs( vec t, mat p)
{
	for (int i = 0; i < p.n_cols; i++)
	{
		coeffs(i) = calcCoeff( t, p.col(i) );
	}
}

void spline::readViaPts( char* fname)
{

	// set T
	mat in;
	in.load( fname );

	setNvp(in.n_rows);	
	setN(in.n_cols-1);

	T = in.col(N);
		
	// actually populate p
	P = in.cols(0,N-1);

}

int spline::getSplineNum( double &t)
{

	t = min( t, T(Nvp-1) );
	// check range on t
	if ( t < T(0) || t > T(Nvp-1) )
	{
		cout << "Spline interpolation: t is out of range" << endl;
		cout << "T(0): " << T(0) << endl;
		cout << "t: " << t << endl;
		cout << "T(N-1): " << T(N-1) << endl;
	}
	int n = 0;
	while ( n < Nvp-1 && ( t < T(n) || t > T(n+1) ) )
	{
		n++;
	}
	return n;
}

double spline::getAQ( int n, int i, double t)
{
	return coeffs(i)(n,0) + coeffs(i)(n,1)*t + coeffs(i)(n,2)*t*t + coeffs(i)(n,3)*t*t*t;
}

/*
int main()
{
	spline sp1( const_cast<char*>("path1.mat") );
	cout << sp1.getNvp() << " " << sp1.getN() << endl;
	sp1.getT().print("T:");
	sp1.getP().print("P:");
	cout << "eval(0,0): " << sp1.eval(1) << endl;
	return 0;
}
*/
