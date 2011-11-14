#ifndef BOX_H
#define BOX_H

using namespace std;
using namespace arma;

class RBbox : public RB
{
	public:
	RBbox(){};
	RBbox( double m, double w, double h, double d)
	{
		RB::setMass(m);
		setW( w);
		setH( h);
		setD( d);
		
		calcCOM();
		calcMCI();
	}
	
	void setW( double w);
	void setH( double h);
	void setD( double d);
	
	double getW();
	double getH();
	double getD();
	
	void calcCOM();
	void calcMCI();
	void draw( vec x0, vec axis);
	
	private:
	double W;
	double H;
	double D;
};

#endif
