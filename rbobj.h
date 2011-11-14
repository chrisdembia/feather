#ifndef RBOBJ_H
#define RBOBJ_H

/*
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
*/
using namespace std;
using namespace arma;

class RBobj : public RB
{
    public:
    RBobj(){};
    RBobj(char* fname)
    {
        objfile= fname;
        readOBJ(objfile);

        calcXC();
        calcYC();
        calcZC();
        calcR();
    };

    void readOBJ(char*);

    int getNumV();
    vector<float> getVX();
    vector<float> getVY();
    vector<float> getVZ();

    float getXC();
    float getYC();
    float getZC();
    double getR();

    int getNumF();
    vector<int> getF1();
    vector<int> getF2();
    vector<int> getF3();

	void normalize();
    void draw( vec pos, vec axis);

    private:
    char* objfile;

    int numV;
    vector<float> VX;
    vector<float> VY;
    vector<float> VZ;

    void calcXC();
    void calcYC();
    void calcZC();
    void calcR();

    float XC;
    float YC;
    float ZC;
    double R;

    int numF;
    vector<int> F1;
    vector<int> F2;
    vector<int> F3;

};

#endif
