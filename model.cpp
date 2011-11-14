#include <fstream>
#include <string>
#include <armadillo>
#include <vector>
#include "spatial.h"
#include "rb.h"

#include "model.h"

using namespace std;
using namespace arma;

void model::setN(int n)
{
	N = n;
	parent.set_size(N);
	pitch.set_size(N);
	Xtree.set_size(N);
	I.set_size(N);
	
	rb.reserve(N);
	
	Xj.set_size(N);
	S.set_size(N);
	vj.set_size(N);
	cj.set_size(N);
}

void model::setParent( int i, int lambda)
{
	parent(i) = lambda;
}

void model::setParent( irowvec pa)
{
	parent = pa;
}

void model::setPitch( int i, double pi)
{
	pitch(i) = pi;
}

void model::setPitch( rowvec p)
{
	pitch = p;
}

void model::setXtree( int i, mat xti)
{
	Xtree(i) = xti;
}

void model::setI( int i, mat ii)
{
	I(i) = ii;
}

void model::setRB( int i, RB* b)
{
	rb[i] = b;
}

void model::setFuture( vec tin, mat future)
{
	t = tin;
	q = future.cols( 0, N-1);
	qd = future.cols( N, 2*N-1);
}

void model::setFuture( vec tin, mat qin, mat qdin)
{
	t = tin;
	q = qin;
	qd = qdin;
}

void model::setTorque( mat tau)
{
	torque = tau;
}

int model::getN()
{
	return N;
}

int model::getParent( int i)
{
	return as_scalar( parent(i) );
}

rowvec model::getPitch()
{
	return pitch;
}

mat model::getXtree( int i)
{
	return Xtree(i);
}

mat model::getI( int i)
{
	return I(i);
}

mat model::getXj( int i)
{
	return Xj(i);
}

mat model::getS( int i)
{
	return S(i);
}

vec model::getVj( int i)
{
	return vj(i);
}

vec model::getCj( int i)
{
	return cj(i);
}

RB* model::getRB( int i)
{
	return rb[i];
}

vec model::getT()
{
	return t;
}

mat model::getQ()
{
	return q;
}

mat model::getQd()
{
	return qd;
}

mat model::getTorque()
{
	return torque;
}

void model::jcalc( int i, double q)
{
	vec s = zeros<vec>(6);

	if ( pitch(i) == 4)
	{
		double rad = getRB(i)->getR();
		vec pos = zeros<vec>(3);
		pos(0) = rad * q;
		Xj(i) = rotz( q ) * xlt( pos );
		s(2) =  1;
		s(3) =  rad * cos(q);
		s(4) = -rad * sin(q);
	}
	else if ( pitch(i) == 0 || pitch(i) == 1 || pitch(i) == 2 || pitch(i) == math::inf() )
	{
		if ( pitch(i) == 0 ) // revolute
		{
			Xj(i) = rotx(q);
			s(0) = 1;
		}
		else if ( pitch(i) == 1 )
		{
			Xj(i) = roty(q);
			s(1) = 1;
		}
		else if ( pitch(i) == 2 )
		{
			Xj(i) = rotz(q);
			s(2) = 1;
		}
		else if ( pitch(i) == math::inf() )
		{
			vec r = zeros<vec>(3);
			r(2) = q;
			Xj(i) = xlt(r);
			s(5) = 1;
		}

	}
	else
	{
		cout << "Joint " << i << " has an invalid jtype." << endl;
	}
	S(i) = s;

}

void model::jcalc( int i, double q, double qd)
{

	vec s = zeros<vec>(6);
	vec VJ = zeros<vec>(6);
	vec CJ = zeros<vec>(6);
	if ( pitch(i) == 4 )
	{

		double rad = getRB(i)->getR();

		vec pos = zeros<vec>(3);
		pos(0) = rad * q;

		Xj(i) = rotz( q ) * xlt( pos );

		s(2) =  1;
		s(3) =  rad * cos(q);
		s(4) = -rad * sin(q);
		VJ(2) =  qd;
		VJ(3) =  rad * cos(q) * qd;
		VJ(4) = -rad * sin(q) * qd;
		CJ(3) = -rad * sin(q) * pow(qd,2);
		CJ(4) = -rad * cos(q) * pow(qd,2);

	}
	else if ( pitch(i) == 0 || pitch(i) == 1 || pitch(i) == 2 || pitch(i) == math::inf() )
	{
		if ( pitch(i) == 0 ) // revolute
		{
			Xj(i) = rotx(q);
			s(0) = 1;
		}
		else if ( pitch(i) == 1 )
		{
			Xj(i) = roty(q);
			s(1) = 1;
		}
		else if ( pitch(i) == 2 )
		{
			Xj(i) = rotz(q);
			s(2) = 1;
		}
		else if ( pitch(i) == math::inf() )
		{
			vec r = zeros<vec>(3);
			r(2) = q;
			Xj(i) = xlt(r);
			s(5) = 1;

		}
		VJ = s * qd;
	}
	else
	{
		cout << "Joint " << i << " has an invalid jtype." << endl;
	}

	S(i) = s;
	vj(i) = VJ;
	cj(i) = CJ;


/* this is for screw transform
		vec r = zeros<vec>(3);
		r(2) = q*pitch(i);
		Xj(i) = rotz(q) * xlt(r);
		s(2) = 1;
		s(5) = pitch(i);
*/
}

void model::makeQFile(char* fname)
{
	int nsteps = q.n_rows;
	int N = q.n_cols;
	ofstream fid( fname );
	fid.precision(10);
	for (int i = 0; i < nsteps; i++)
	{
		fid << t(i);
		for (int j = 0; j < N; j++)
		{
			fid << " " << q(i,j);
		}
		for (int j = 0; j < N; j++)
		{
			fid << " " << qd(i,j);
		}
		fid << endl;
	}
	fid.close();
}


void model::printModel()
{
	cout << "N: " << N << endl;
	cout << endl;
	
	cout << "parent:" << endl;
	cout << parent << endl;
	cout << endl;
	
	cout << "pitch:" << endl;
	cout << pitch << endl;
	
	for (int i = 0; i < N; i++)
	{
		cout << "Xtree(" << i << "):" << endl;
		cout << Xtree(i) << endl;
	}
	cout << endl;
	
	for (int i = 0; i < N; i++)
	{
		cout << "I(" << i << "):" << endl;
		cout << I(i) << endl;
	}
	cout << endl;
	
	cout << "END OF DISPLAYMODEL -- END OF DISPLAYMODEL" << endl;
}

void model::saveFuture( string fname)
{
	t.save( fname + "_t.mat", raw_ascii);
	q.save( fname + "_q.mat", raw_ascii);
	qd.save( fname + "_qd.mat", raw_ascii);
}

void model::loadFuture( string fname)
{
	t.load( fname + "_t.mat", raw_ascii);
	q.load( fname + "_q.mat", raw_ascii);
	q.print();
	qd.load( fname + "_qd.mat", raw_ascii);
	qd.print();
}

void model:: saveTorque( string fname )
{
	torque.save( fname + "_tau.mat", raw_ascii);
}
