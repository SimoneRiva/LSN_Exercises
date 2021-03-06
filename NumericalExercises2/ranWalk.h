#ifndef __RANWALK__
#define __RANWALK__
#include <cmath>
#include <iostream>

using namespace std;

class Vector3D {

protected:

	double * _x = new double [3];

public:

	Vector3D() {_x[0]=0.; _x[1]=0.; _x[2]=0.;};
	Vector3D(double x, double y, double z) {_x[0]=x; _x[1]=y; _x[2]=z;};
	~Vector3D();

//	Vector3D& operator= (const Vector3D& v);

	void SetCartesianCoords(double x, double y, double z) {_x[0]=x; _x[1]=y; _x[2]=z;};
	void SetSphericalCoords(double r, double theta, double phi) {_x[0]=r*sin(theta)*cos(phi); _x[1]=r*sin(theta)*sin(phi); _x[2]=r*cos(theta);};
	double X() {return _x[0];}
	double Y() {return _x[1];}
	double Z() {return _x[2];}
	double theta() {return acos(_x[2]/R());}
	double phi() {return atan(_x[1]/_x[0]);}
	double R() {return sqrt(_x[0]*_x[0]+_x[1]*_x[1]+_x[2]*_x[2]);}
	double squareR() {return (_x[0]*_x[0]+_x[1]*_x[1]+_x[2]*_x[2]);}

  void increase(int d, double value); //Increases d-th component of the Vector3D by value

	void discRWStep(double a, double r, double s); //Modifies the Vector3D into the following step of a random walk on a cubic lattice of constant a (r and s must be provided  uniformly distributed in [0,3) and [0,2) respectively).

	void contRWStep(double a, double theta, double phi); //Modifies the Vector3D into the following step (length a) of a random walk in continuum space (theta and phi must be provided such that they sample the solid angle uniformly).

};

#endif
