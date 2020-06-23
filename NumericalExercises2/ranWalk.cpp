#include <cmath>
#include <iostream>
#include "ranWalk.h"

using namespace std;

/*
Vector3D& Vector3D::operator=(const Vector3D& v) {
	_x[0]=v.X();
	_x[1]=v.Y();
	_x[2]=v.Z();
	return *this;
}
*/

void Vector3D::increase(int d, double value){
	if(d>2) cerr << "Error: Vector3d only has 3 components!" << endl;
	_x[d] += value;
}

void Vector3D::discRWStep(double a, double r, double s){
	this->increase((int)r, a*(2*(int)s-1));
}

void Vector3D::contRWStep(double a, double theta, double phi){
	this->increase(0, a*sin(theta)*cos(phi));
	this->increase(1, a*sin(theta)*sin(phi));
	this->increase(2, a*cos(theta));
}
