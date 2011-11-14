#ifndef MODEL_H
#define MODEL_H
/*
#include <armadillo>
#include <vector>
#include "spatial.h"

*/
#include "rb.h"
//#include "rbbox.h"
//#include "rbobj.h"

using namespace std;
using namespace arma;

// 3D, not planar
class model
{
	public:
	model(){};
	model(int n, irowvec pa, field<mat> xt, rowvec pi, field<mat> i)
	{
		setN(n);
		setParent(pa);
		Xtree = xt;
		setPitch(pi);
		I = i;
	}
	
	// setters
	void setN(int);
	void setParent( int i, int lambda);
	void setParent( irowvec pa);
	void setPitch( int i, double pi);
	void setPitch( rowvec p);
	void setXtree( int i, mat xti);
	void setI( int i, mat ii);
	void setRB( int i, RB* b);
	void setFuture( vec tin, mat future);
	void setFuture( vec tin, mat q, mat qd);
	void setTorque( mat tau);
	
	// getters
	int getN();
	int getParent( int i);
	rowvec getPitch();
	mat getXtree( int i);
	mat getI( int i);
	mat getXj( int i);
	mat getS( int i);
	vec getVj( int i);
	vec getCj( int i);
	RB* getRB( int i);
	
	vec getT();
	mat getQ();
	mat getQd();
	mat getTorque();
	
	void jcalc( int i, double q);
	void jcalc( int i, double q, double qd);	
	void makeQFile(char* fname);
	void printModel();
	void saveFuture( string fname);
	void loadFuture( string fname);
	void saveTorque( string fname);
	
	private:
	int N;
	irowvec parent;
	rowvec pitch;
	field<mat> Xtree;
	field<mat> I;
	vector<RB*> rb;

	field<mat> Xj;
	field<mat> S;
	field<vec> vj;
	field<vec> cj;
	
	vec t;
	mat q;
	mat qd;
	mat torque;
};


#endif
