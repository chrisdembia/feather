#include <armadillo>
#include "model.h"
#include "spatial.h"

#include "featherstone.h"

using namespace std;
using namespace arma;

vec FD( model &M, vec q, vec qd, vec tau)
{
	const int N = M.getN();
	
	vec qdd(N);
	
	// set acceleration of base
	vec a_base = zeros<vec>(6);
	a_base(4) = -9.81;
	
	field<vec> S(N);
	field<mat> Xup(N);
	field<vec> v(N);
	field<vec> c(N);
	field<mat> IA(N);
	field<vec> pA(N);
	
	vec vj(6);

	for (int i = 0; i < N; i++)
	{
		M.jcalc( i, q(i), qd(i) );
		S(i) = M.getS(i);
		Xup(i) = M.getXj( i) * M.getXtree(i);
		
		vec vj = M.getVj(i);
		vec cj = M.getCj(i);

		if ( M.getParent(i) == 0 )
		{
			v(i) = vj;
			c(i) = cj;
		}
		else
		{
			v(i) = Xup(i) * v(M.getParent(i)-1) + vj;
			c(i) = crm( v(i) ) * vj + cj;
		}

		IA(i) = M.getI(i);
		pA(i) = crf( v(i) ) * M.getI(i) * v(i);
		
		// add in external force to pA at this point
	} // end of forward iteration

	field<vec> U(N);
	vec d(N); // THIS WILL NOT WORK FOR GENERAL JOINTS, THAT ARE NOT 6X1 S
	vec u(N);
	
	mat Ia(6,6);
	vec pa(6);

	for (int i = N-1; i >= 0; i--)
	{
		U(i) = IA(i) * S(i);
		d(i) = dot( S(i), U(i));
		u(i) = tau(i) - dot( S(i), pA(i));
		if ( M.getParent(i) != 0 )
		{
			Ia = IA(i) - U(i) / d(i) * trans(U(i));
			pa = pA(i) + Ia * c(i) + U(i) * u(i)/d(i);
			IA(M.getParent(i)-1) = IA(M.getParent(i)-1) + trans(Xup(i)) * Ia * Xup(i);
			pA(M.getParent(i)-1) = pA(M.getParent(i)-1) + trans(Xup(i)) * pa;
		}
	}
	
	field<vec> a(N);

	for (int i = 0; i < N; i++)
	{
		if (M.getParent(i) == 0)
		{
			a(i) = Xup(i) * -a_base + c(i);
		}
		else
		{
			a(i) = Xup(i) * a(M.getParent(i)-1) + c(i);
		}
		qdd(i) = (u(i) - dot( U(i), a(i)) ) / d(i);
		a(i) = a(i) + S(i) * qdd(i);
	}
	return qdd;
}

mat bodypos( model &M, int b, vec q)
{
	mat X = eye<mat>(6,6);
	mat XJ(6,6);
	
	while ( b >= 0 )
	{
		M.jcalc( b, q(b) );
		X = X * M.getXj(b) * M.getXtree(b);
		b = M.getParent(b)-1;
	}
	
	return X;
}

vec xformpos( mat X)
{
	vec pos(3);

	mat rcrosst = trans( X( span(0,2), span(0,2) ) ) * X( span(3,5), span(0,2) );
	
	pos(0) = rcrosst(1,2);
	pos(1) = rcrosst(2,0);
	pos(2) = rcrosst(0,1);
	
	return pos;
}

vec xformrot ( mat X)
{
	vec v(3);
	
	mat E = X.submat(0,0,2,2);
	
	double RT_RESOL = 4.5e-9;
	double SMALL = 1e-4;
	
	v(0) = E(1,2) - E(2,1);
	v(1) = E(2,0) - E(0,2);
	v(2) = E(0,1) - E(1,0);
	
	double S = dot( v, v);
	
	double tr1 = trace(E) - 1.0;
	
	if ( S < RT_RESOL && tr1 > 0 )
	{
		S = ( S + 24 ) / 48;
		v = v*S;
	}
	else if ( S < SMALL && tr1 < 0 )
	{
		S = atan2( sqrt(S), tr1 );
		
		double d;
		
		if ( E(0,0) > E(1,1) && E(0,0) > E(2,2) )
		{
			d = 2 * E(0,0) - tr1;
			v(0) = S * sqrt( d / (2 - tr1) );
			
			if ( E(2,1) > E(1,2) )
			{
				v(0) = -v(0);
			}
			
			v(1) = ( E(0,1) + E(1,0) ) * v(0) / d;
			v(2) = ( E(0,2) + E(2,0) ) * v(0) / d;
		}
		else if ( E(1,1) > E(2,2) )
		{
			d = 2 * E(1,1) - tr1;
			v(1) = S * sqrt( d / (2 - tr1) );
			
			if ( E(0,2) > E(2,0) )
			{
				v(1) = -v(1);
			}
			
			v(0) = ( E(0,1) + E(1,0) ) * v(1) / d;
			v(2) = ( E(1,2) + E(2,1) ) * v(1) / d;
		}
		else
		{
			d = 2 * E(2,2) - tr1;
			v(2) = S * sqrt( d / (2 - tr1) );
			
			if ( E(1,0) > E(0,1) )
			{
				v(2) = -v(2);
			}
			
			v(0) = ( E(0,2) + E(2,0) ) * v(2) / d;
			v(1) = ( E(1,2) + E(2,1) ) * v(2) / d;			
		}
	}
	else
	{
		S = sqrt(S);
		S = atan2( S, tr1) / S;
		
		v = v * S;
	}
	
	return v;
}

vec rel2abs( model &M, int b, vec q, vec relpos)
{
	mat X = bodypos( M, b, q);

	vec abspos = xformpos( X) + inv( X.submat(0,0,2,2) ) * relpos;

	return abspos;
}

vec rel2abschainp( model &M, int b, vec q, vec relpos)
{
	vec abspos = zeros<vec>(3);
	vec e = zeros<vec>(3);
	e(0) = 1; /// should be length of member
	for (int i = 0; i < b; i++)
	{
		abspos += inv( rz( sum( q.rows(0,i) ) ) )*e;
	}
	return abspos;
}

vec abs2rel( model &M, int b, vec q, vec abspos)
{
	mat X = bodypos( M, b, q);
	vec relpos = X( span(0,2), span(0,2) ) * abspos;
	relpos(0) += X(5,1);
	relpos(1) += X(3,2);
	relpos(2) += X(4,0);

	return relpos;
}

mat bodyjac( model &M, int b, vec q)
{
	int N = M.getN();
	rowvec e = zeros<rowvec>(N);

	while ( b >= 0 )
	{
		e(b) = 1;
		b = M.getParent(b)-1;
	}
	
	mat J = zeros<mat>(6,N);

	field<mat> Xa(N);	
	for (int i = 0; i < N; i++)
	{
		if ( e(i) )
		{
			M.jcalc( i, q(i) );
			Xa(i) = M.getXj(i) * M.getXtree(i);

			if ( M.getParent(i) != 0 )
			{
				Xa(i) = Xa(i) * Xa( M.getParent(i)-1 );
			}
			J.col(i) = solve( Xa(i) , M.getS(i) );
		}
	}
	return J;
}

field<vec> IKpos( model &M, int b, mat Xd, vec q0)
{
	
	double lambda = 0.0001;

	field<vec> q(1);
	q(0) = q0;
	field<vec> qtemp;
	vec dpos = ones<vec>(6);
	vec dq(6);
	mat X(6,6);
	mat J;
	int i = 0;

	while ( norm( dpos , 2) > 1e-5 )
	{
		X = bodypos( M, b, q(i) );
		J = bodyjac( M, b, q(i) );
		J = X * J;

//		dpos = XtoV( solve( Xd , X ) );
		dpos = XtoV( Xd * inv(X) );
		dpos.print("dpos:");
//		dq = solve( trans(J) * J, trans(J) * dpos);
//		dq = trans(J) * solve( ( J * trans(J) ) , dpos );
//		dq = solve( trans(J) * J + pow(lambda,2)*eye(M.getN(),M.getN()) , trans(J) * dpos );
		dq = trans(J) * solve( J * trans(J) + pow(lambda,2) * eye(6,6), dpos );

		i++;

		qtemp = q;

		q.set_size(i+1);

		q.rows(0,i-1) = qtemp;

		q(i) = q(i-1) + dq;
	}
	return q;
}





