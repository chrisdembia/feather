#include <iostream>
#include <cstdlib>
#include <armadillo>
#include <vector>
#include <string>
#include "spatial.h"
#include "rb.h"
#include "rbbox.h"
#include "rbobj.h"
#include "rbnull.h"
#include "model.h"
#include "spline.h"
#include "prophet.h"
#include "featherstone.h"
#include "play.h"

using namespace std;
using namespace arma;

void chain( model &M, int N);
void chainrot( model &M, int nlinks);
void chain3DOF( model &M);
void chain3DOF( model &M, int nlinks);
void chain3DOF2( model &M, int nlinks);

int main()
{

	int doODE = 1;	
	string matname = "accordion0430";
	
	wall_clock runtime;
	runtime.tic();

	cout << "Creating model" << endl;

	model p;

//	int N = 10;
//	chain( p, N);
	
	int nlinks = 10;
	chainrot( p, nlinks);
	int N = p.getN();
	
//	chain3DOF( p);
//	int N = p.getN();

//	int nlinks = 5;
//	chain3DOF( p, nlinks);
//	chain3DOF2( p, nlinks);
//	int N = p.getN();

	//p.displayModel();

	int FPS = 30; // frames per second
	int T = 20; // duration of animation
	vec tspan = linspace<vec>( 0.0, double(T), FPS*T); // the number of frames we want

	// PATH PLANNING
	spline sp1( const_cast<char*>("path4accordion.mat") );

	// INITIAL CONDITIONS
	vec q0 = zeros<vec>(N);
//	q0 = sp1.eval(0);
	vec qd0 = zeros<vec>(N);
//qd0(2) = 10.0f;
//	q0(0) = -math::pi()/2;	
//q0(0) = math::pi()/4;
//qd0(0) = 0.5f;
//qd0(1) = 0.0f;
//qd0(2) = 1.0f;
	vec ICs(2*N);
	ICs.rows(0,N-1) = q0;
	ICs.rows(N,2*N-1) = qd0;

	// CONTROL
	vec Ks = zeros<vec>(2*N);
//	vec Ks = 1*ones<vec>(2*N);
	Ks.rows(0,N-1) = linspace<vec>(10000,5000, N);
	Ks.rows(N,2*N-1) = linspace<vec>(1000,500, N);
	
	// SOLVE
	if ( doODE )
	{
		cout << "Simulating..." << endl;
		mat Future = prophet( p, tspan, ICs, Ks, sp1);
		cout << "Storing simulation in model" << endl;
		p.setFuture( tspan, Future);
		cout << "Storing simulation in ARMA mat files" << endl;
		p.saveFuture( matname);
	}
	else
	{
		cout << "Loading simulation from ARMA mat files" << endl;
		p.loadFuture( matname);
	}
	
	cout << "Storing simulation in output file" << endl;
	char* tablet = const_cast<char*>("tablet.dat");
	p.makeQFile(tablet);
	
	// STORE TORQUE
	vec qtemp(N);
	vec qktemp(N);
	vec qdtemp(N);
	vec tautemp(N);
	mat torques(FPS*T, N);
	mat QK(FPS*T,N);
		
	for (int i = 0; i < FPS*T; i++ )
	{	
		qtemp = vec( trans( p.getQ().row(i) ) );
		qktemp = sp1.eval( p.getT()(i) );
		qdtemp = vec( trans ( p.getQd().row(i) ) );
		
		tautemp = -(Ks.rows( 0, N-1) % ( qtemp - qktemp ) ) - (Ks.rows( N, 2*N-1) % qdtemp );
		torques.row(i) =  rowvec( trans( tautemp ) );
		
		QK.row(i) = rowvec( trans( qktemp ) );
	}
	p.setTorque( torques);
	p.saveTorque( matname);
	QK.save( matname + "_qk.mat", raw_ascii);
	
	// PLAY ANIMATION
	int saveImgs = 0;
	double rauntime = runtime.toc();
	cout << "Simulation took " << runtime.toc() << " s to prepare." << endl;
	int k = play( p, sp1, saveImgs);
	cout << "Duration of simulation: " << runtime.toc() - rauntime << " s" << endl;

	return 0;
}

void chain( model &M, int N)
{
	cout << "Creating a chain model" << endl;
	// make a chain model

	double m = 1.0; // mass of links
	double L = 1.0; // length of each chain link
	
	// define N, parent, pitch, Xtree, I
	M.setN(N);
	
	// deifne size
	if (N == 1)
	{
		M.setParent( zeros<irowvec>(1) );
	}
	else
	{
		M.setParent( linspace<irowvec>( 0, N-1, N) );
	}
	
	// define vector (3D pos) argument for xtree
	vec r = zeros<vec>(3);
	M.setXtree( 0, xlt(r) );

	// cdefine center of mass 3D pos vector
	vec com = zeros<vec>(3);
	com(0) = L/2.0;

	// define inertia
	vec inertia = m*pow(L,2)/12.0 * ones<vec>(3);
	inertia(0) = 0.01; // no thickness in this direction
	M.setI(0, mcI( m, com, diagmat(inertia)) );
		
	// RIGID BODY
	RBbox* box1 = new RBbox( m, L, L/5, L/5);

/*	char* torus = const_cast<char*>("torus.obj");
	RBobj* box1 = new RBobj( torus );
	box1->normalize();
*/
	M.setRB( 0, box1);

	r(0) = L; // for the rest of the links
	// fill in the fields
	for (int i = 1; i < N; i++)
	{
		M.setXtree( i, xlt(r) );
		M.setI( i, mcI( m, com, diagmat(inertia)) );
		M.setRB( i, box1);
	}
	
	// define pitch
	M.setPitch( 2*ones<rowvec>(N) );

}

void chainrot( model &M, int nlinks)
{
	cout << "Creating a chainrot model" << endl;
	// make a chain model

	int N = nlinks + 1;
	
	double m = 1.0; // mass of links
	double L = 1.0; // length of each chain link
	
	// define N, parent, pitch, Xtree, I
	M.setN(N);
	
	// deifne size
	M.setParent( linspace<irowvec>( 0, N-1, N) );
	
	// define vector (3D pos) argument for xtree
	vec r = zeros<vec>(3);
	// RIGID BODY
	RBnull* null1 = new RBnull();
	RBbox* box1 = new RBbox( m, L, L/6, L/6);
	
	// cdefine center of mass 3D pos vector
	vec com = zeros<vec>(3);
	com(0) = L/2.0;

	// define inertia
	vec inertia = m*pow(L,2)/12.0 * ones<vec>(3);
	inertia(0) = 0.01; // no thickness in this direction
	
	M.setXtree( 0, xlt(r) );
	M.setI(0, mcI( 0, zeros<vec>(3), diagmat(zeros<vec>(3) ) ) );
	M.setRB( 0, null1);
		
	M.setXtree( 1, xlt(r) );
	M.setI(1, mcI( m, com, diagmat(inertia)) );
	M.setRB( 1, box1);

	r(0) = L; // for the rest of the links
	// fill in the fields
	for (int i = 2; i < N; i++)
	{
		M.setXtree( i, xlt(r) );
		M.setI( i, mcI( m, com, diagmat(inertia)) );
		M.setRB( i, box1);
	}
	
	// define pitch
	rowvec pitch = 2*ones<rowvec>(N);
	pitch(0) = 1;
	M.setPitch( pitch );

}

void chain3DOF( model &M)
{
	cout << "Creating a chain3DOF model" << endl;
	// make a chain model

	double m = 0.0; // mass of links
	double L = 1.0; // length of each chain link
	
	int N = 3;
	
	// define N, parent, pitch, Xtree, I
	M.setN(N);
	
	// deifne size
	M.setParent( linspace<irowvec>( 0, N-1, N) );
	
	// define vector (3D pos) argument for xtree
	M.setXtree(0,  xlt( zeros<vec>(3) ) );
	M.setXtree(1, rotx( -math::pi()/2 ) );
	M.setXtree(2, roty( math::pi()/2 ) );

	// cdefine center of mass 3D pos vector
	vec com = zeros<vec>(3);

	vec inertia = zeros<vec>(3);
	
	M.setI(0, mcI( m, com, diagmat(inertia)) );
	M.setI(1, mcI( m, com, diagmat(inertia)) );

	com(0) = L/2.0;

	m = 1.0;
	
	// define inertia
	inertia = m*pow(L,2)/12.0 * ones<vec>(3);
	inertia(0) = 0.01; // no thickness in this direction
	M.setI(2, mcI( m, com, diagmat(inertia)) );
		
	// RIGID BODY
	RBbox* box1 = new RBbox( m, L, L/5, L/5);
	RBnull* null1 = new RBnull();
	
	M.setRB( 0, null1);
	M.setRB( 1, null1);
	M.setRB( 2, box1);
	
	// define pitch
	M.setPitch( 2*ones<rowvec>(N) );

/*
	M.setXtree(0, xlt(zeros<vec>(3)));
	M.setXtree(1, rotx( math::pi()/2) );
	M.setXtree(2, roty( math::pi()/2) );
*/
	

}

void chain3DOF( model &M, int nlinks)
{
	cout << "Creating a chain3DOF model" << endl;
	// make a chain model

	double m = 0.0; // mass of links
	double L = 1.0; // length of each chain link
	
	int N = 3*nlinks;
	
	// define N, parent, pitch, Xtree, I
	M.setN(N);
	
	// deifne size
	M.setParent( linspace<irowvec>( 0, N-1, N) );
	
	// define vector (3D pos) argument for xtree
	vec com;
	vec inertia;

	RBbox* box1 = new RBbox( m, L, L/5, L/5);
	RBnull* null1 = new RBnull();

	M.setXtree(0,  xlt( zeros<vec>(3) ) );
	vec r(3);
	r(0) = L;
	for (int i = 0; i < N-2; i+=3)
	{
		if ( i != 0 )
		{
				M.setXtree(0+i,  xlt( r ) );
		}
		M.setXtree(1+i, rotx( -math::pi()/2 ) );
		M.setXtree(2+i, roty( math::pi()/2 ) );
		//
		m = 0.0;
		com = zeros<vec>(3);
		inertia = zeros<vec>(3);
		M.setI(0+i, mcI( m, com, diagmat(inertia)) );
		M.setI(1+i, mcI( m, com, diagmat(inertia)) );
		
		com(0) = L/2.0;	m = 1.0;
		inertia = m*pow(L,2)/12.0 * ones<vec>(3); inertia(0) = 0.01;
		M.setI(2+i, mcI( m, com, diagmat(inertia)) );	
		//
		M.setRB( 0+i, null1);
		M.setRB( 1+i, null1);
		M.setRB( 2+i, box1);
			
	}
	// define pitch
	M.setPitch( 2*ones<rowvec>(N) );

}

void chain3DOF2( model &M, int nlinks)
{
	cout << "Creating a chain3DOF2 model" << endl;
	// make a chain model

	double m = 0.0; // mass of links
	double L = 1.0; // length of each chain link
	
	int N = 3*nlinks;
	
	// define N, parent, pitch, Xtree, I
	M.setN(N);
	
	// deifne size
	M.setParent( linspace<irowvec>( 0, N-1, N) );
	
	// define vector (3D pos) argument for xtree
	vec com;
	vec inertia;

	RBbox* box1 = new RBbox( m, L, L/5, L/5);
	RBnull* null1 = new RBnull();

	M.setXtree(0,  xlt( zeros<vec>(3) ) );
	vec r(3);
	r(0) = L;
	
	for (int i = 0; i < N-2; i+=3)
	{
		if ( i != 0 )
		{
				M.setXtree(0+i,  xlt( r ) );
		}
		M.setXtree(1+i, xlt( zeros<vec>(3) ) );
		M.setXtree(2+i, xlt( zeros<vec>(3) ) );
		//
		m = 0.0;
		com = zeros<vec>(3);
		inertia = zeros<vec>(3);
		M.setI(0+i, mcI( m, com, diagmat(inertia)) );
		M.setI(1+i, mcI( m, com, diagmat(inertia)) );
		
		com(0) = L/2.0;	m = 1.0;
		inertia = m*pow(L,2)/12.0 * ones<vec>(3); inertia(0) = 0.01;
		M.setI(2+i, mcI( m, com, diagmat(inertia)) );	
		//
		M.setRB( 0+i, null1);
		M.setRB( 1+i, null1);
		M.setRB( 2+i, box1);
		
		M.setPitch(i+0, 0);			
		M.setPitch(i+1, 1);			
		M.setPitch(i+2, 2);			
	}
	// define pitch

}
/*
void bike( model &M)
{
  cout << "Creating a bike model" << endl;
  // use the benchmark parameters
  int N = 4;
  M.setN(4);
  M.setParent( linspace<irowvec>( 0, N-1,N) );
}
*/

/*	// bodyJac 
	vec q02 = zeros<vec>(N);
	
	// IKpos 
	vec poss = zeros<vec>(3);
	poss(1) = N-6;
//	poss(1) = N-1;
//	poss(1) = 0.7071;
	mat Xd =  rotz( math::pi() / 4 ) * xlt( poss ); 
	
//	q02(0) = 0.01; //math::pi()/2;
//	q02(1) = 120*math::pi()/180;
//	q02(2) = 0;
	field<vec> qout = IKpos( p, N-5, Xd, q02);
	poss(1) = 1;
	poss(0) = 1;
	Xd = xlt( poss );
	qout = IKpos( p, N-1, Xd, qout(qout.n_elem-1) );
	
	cout << "qout: " << qout << endl;
	cout << "number of iterations: " << qout.n_elem << endl;
	int l = play( p, qout);
*/	

/*	field<vec> qk(3);
	qk(0) = sp1.eval( 0 );
	qk(1) = sp1.eval( 1 );
	cout << qk(1) << endl;
	qk(2) = sp1.eval( 2 );
	int l = play( p, qk, 0);
*/
