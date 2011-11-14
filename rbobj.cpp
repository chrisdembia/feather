#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "rb.h"

#include "rbobj.h"

using namespace std;
using namespace arma;

void RBobj::readOBJ(char* objfile)
{
    string line;
	numV = 0;
	numF = 0;

	ifstream objh( objfile );
	if (objh.is_open())
	{
		while(!objh.eof())
		{
			getline(objh,line);
			if (line.compare(0,1,"v") == 0)
			{
				numV++;
			}

			if (line.compare(0,1,"f") == 0)
			{
				numF++;
			}
		}
	}
	objh.close();

    VX.reserve(numV);
    VY.reserve(numV);
    VZ.reserve(numV);
    F1.reserve(numF);
    F2.reserve(numF);
    F3.reserve(numF);

    string garbage1;
	string temp1;
	string temp2;
	string temp3;

	int vidx = 0;
	int fidx = 0;
	ifstream objh2( objfile );
	if (objh2.is_open())
	{
		while(!objh2.eof())
		{
			objh2 >> line;
			if (line.compare(0,1,"#") == 0)
			{
				getline(objh2,line);
			}
			else
			{
				objh2 >> temp1 >> temp2 >> temp3;
				if (line.compare(0,1,"v") == 0)
				{
                    VX.push_back(atof(temp1.c_str()));
                    VY.push_back(atof(temp2.c_str()));
                    VZ.push_back(atof(temp3.c_str()));

					//cout << "v " << vidx << " " << vx[vidx] << " " << vy[vidx] << " " << vz[vidx] << endl;
					vidx++;
				}
				else if (line.compare(0,1,"f") == 0)
				{
                    F1.push_back(atof(temp1.c_str()));
                    F2.push_back(atof(temp2.c_str()));
                    F3.push_back(atof(temp3.c_str()));

					//cout << "f " << fidx << " " << F1[fidx] << " " << F2[fidx] << " " << F3[fidx] << endl;
					fidx++;
				}
			}
		}
	}

	objh2.close();
}

int RBobj::getNumV()
{
    return numV;
}

vector<float> RBobj::getVX()
{
    return VX;
}

vector<float> RBobj::getVY()
{
    return VY;
}

vector<float> RBobj::getVZ()
{
    return VZ;
}

float RBobj::getXC()
{
    return XC;
}

float RBobj::getYC()
{
    return YC;
}

float RBobj::getZC()
{
    return ZC;
}

double RBobj::getR()
{
	return R;
}

int RBobj::getNumF()
{
    return numF;
}

vector<int> RBobj::getF1()
{
    return F1;
}

vector<int> RBobj::getF2()
{
    return F2;
}

vector<int> RBobj::getF3()
{
    return F3;
}

void RBobj::calcXC()
{
    float sumX = 0;
    for (int i = 0; i < numV; i++)
    {
        sumX += VX[i];
    }
    XC = sumX/float(numV);
}

void RBobj::calcYC()
{
    float sumY = 0;
    for (int i = 0; i < numV; i++)
    {
        sumY += VY[i];
    }
    YC = sumY/float(numV);
}

void RBobj::calcZC()
{
    float sumZ = 0;
    for (int i = 0; i < numV; i++)
    {
        sumZ += VZ[i];
    }
    ZC = sumZ/float(numV);
}

void RBobj::calcR()
{
	float sumR = 0;
	for (int i = 0; i < numV; i++)
	{
		sumR += sqrt( pow( (VX[i] - XC),2) + pow( (VY[i] - YC),2) + pow( (VZ[i] - ZC) ,2) );
	}
	
	R = sumR/(float)numV;
//	R = 13.5;
}

void RBobj::normalize()
{
	for (int i = 0; i < numV; i++)
	{
		VX[i] = (VX[i] - XC) / (2*R); // + 0.5;
		VY[i] = (VY[i] - YC) / (2*R);
		VZ[i] = (VZ[i] - ZC) / (2*R);
	}
	vector<double> znew;
	znew.reserve(numV);
	for (int i = 0; i < numV; i++)
	{
		znew[i] = -VY[i];
		VY[i] = VZ[i];
		VZ[i] = znew[i];
	}
	znew.clear();
	R = 0.5;
}

void RBobj::draw( vec pos, vec axis)
{
	glTranslatef( pos(0), pos(1), pos(2));
	glRotatef( 180/math::pi() * sqrt( dot( axis, axis) ) , axis(0), axis(1), axis(2) );
    for (int i = 0; i < numF; i++ )
	{
		glBegin(GL_TRIANGLES);
            glColor3f((float)i/(float)numF,0.0f,(float)(numF-i)/(float)numF);
/*
            glVertex3f(VX[F1[i]-1]-XC,VY[F1[i]-1]-YC,VZ[F1[i]-1]-ZC);
            glVertex3f(VX[F2[i]-1]-XC,VY[F2[i]-1]-YC,VZ[F2[i]-1]-ZC);
            glVertex3f(VX[F3[i]-1]-XC,VY[F3[i]-1]-YC,VZ[F3[i]-1]-ZC);
*/

            glVertex3f(VX[F1[i]-1],VY[F1[i]-1],VZ[F1[i]-1]);
            glVertex3f(VX[F2[i]-1],VY[F2[i]-1],VZ[F2[i]-1]);
            glVertex3f(VX[F3[i]-1],VY[F3[i]-1],VZ[F3[i]-1]);
		glEnd();
	}
}
