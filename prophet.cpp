#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <armadillo>
#include "model.h"
#include "featherstone.h"
#include "spline.h"

#include "prophet.h"

using namespace std;
using namespace arma;

struct PARAM
{
	model* mptr;
	vec K;
	spline* sptr;
};
     
int crystalball(double t, const double CQ[], double Cdots[], void *params)
{
	// read in some data, yo
	PARAM p = *(PARAM *)params;
	model* Mptr = p.mptr;
	vec Ks = p.K;
	spline* Sptr = p.sptr;
	
	int N = Ks.n_elem/2;

	// prepare velocity and acceleration
    vec Q = CtoARMA( const_cast<double*>(CQ), 2*N);
    vec q = Q.rows( 0, N-1);
    vec qd = Q.rows( N, 2*N-1);
    
    // calculate control torques
//	vec jtorques = -(Ks.rows( 0, N-1) % q) - (Ks.rows( N, 2*N-1) % qd);	

	vec jtorques = -(Ks.rows( 0, N-1) % (q - Sptr->eval(t) ) ) - (Ks.rows( N, 2*N-1) % qd);
    // get derivatives
	vec qdd = FD( *Mptr, q, qd, jtorques);

	// put derivatives next to each other
	vec dots(2*N);
	dots.rows(0, N-1) = qd;
	dots.rows( N, 2*N-1) = qdd;

	ARMAtoC( dots, Cdots);
	    
    return GSL_SUCCESS;
}
     
mat prophet(model &M, vec tspan, vec ICs, vec Ks, spline &sp1)
{

	int nsteps = tspan.n_elem;	
	
	// prepare ICs
	int ndim = ICs.n_elem; // number of equations to solve
	double CQ[ndim]; // initial conditions
	
	ARMAtoC(ICs, CQ);
	
	// prepare parameters for crystalball
	PARAM p1;
	p1.mptr = &M;
	p1.K = Ks;
	p1.sptr = &sp1;
	
	// PREPARE OUTPUT
	mat Future(nsteps,ndim);

	// GSL FUNKS
	cout << "Starting GSL integrator" << endl;
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, ndim);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-4, 1e-4);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (ndim);
	gsl_odeiv_system sys = {crystalball, NULL, ndim, &p1};
	double t = 0.0;
	double ti;
	double h = 1e-6;
	int maxIt = 1000;
	int it;
	
    for (int i = 0; i < nsteps; i++)
    {
    	ti = tspan(i);
    	
    	if ( fmod(ti,1) < 0.032 )
    	{
    		cout << "Integrating up to t = " << ti << endl;
    	}
    	
    	it = 0;
		while ( t < ti )
		{
			it++;
			int status = gsl_odeiv_evolve_apply (e, c, s,
		                                            &sys, 
		                                            &t, ti,
		                                            &h, CQ);
		                                            
		    if ( it >= maxIt )
		    {
	    		cout << "PROPHET: This is taking too long: " << it << " iterations" << endl;
	    		cout << "ti = " << ti << ", t = " << t << endl;
	    		break;
	    	}
		    
			if (status != GSL_SUCCESS)
			{
				cout << "PROPHET: Integration Error" << endl;
				break;
			}
		}
		
		for (int j = 0; j < ndim; j++)
		{
			Future(i,j) = CQ[j];
		}
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return Future;
}

void ARMAtoC( vec v, double y[])
{
	for (int i = 0; i < v.n_elem; i++)
	{
		y[i] = v(i);
	}
}
	
vec CtoARMA( double v[], int N)
{
	vec y(N);
	for (int i = 0; i < N; i++)
	{
		y(i) = v[i];
	}
	return y;
}
/*
	// process output file
	Out(0).set_size(nsteps,1);
	Out(1).set_size(nsteps,ndim);
	
	ifstream fid2( tablet );
	for (int tidx = 0; tidx < nsteps; tidx++)
	{
		fid2 >> Out(0)(tidx);
		for (int dimidx = 0; dimidx < ndim; dimidx++)
		{
			fid2 >> Out(1)(tidx,dimidx);

		}

	}

	fid2.close();

*/

