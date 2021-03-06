#ifndef _ex05_h_
#define _ex05_h_
#include <cmath>

Random rnd;

const int nstep=2e6, equilstep=1e4, nblocks=200, blocksize=nstep/nblocks;
int step, accepted; //To evaluate the acceptance rate of Metropolis
double l, sigma; //Parameters for the transition probability T(x|y)
double x[3];

//Quantities for data blocking
double r[nblocks], ave_r[nblocks], ave2_r[nblocks], err_r[nblocks];
double sum;

ofstream out_points, out_radii;

//Classes
class BaseFunction3D {
	public:
		virtual double Eval(double * x) const=0;
};

class square_1s : public BaseFunction3D {
	public:
		virtual double Eval(double * x) const {return exp(-2*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/M_PI; }
};

class square_2p : public BaseFunction3D {
	public:
		virtual double Eval(double * x) const {
			double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
			return 2*exp(-r)*x[2]*x[2]/(64*M_PI);
		}
};


class MetroStep {
	public:
		virtual void Step(BaseFunction3D * p) const =0;
};

class MetroStepUnif : public MetroStep {
	public:
		virtual void Step(BaseFunction3D * p) const;
};

class MetroStepGauss : public MetroStep {
	public:
		virtual void Step(BaseFunction3D * p) const;
};

//Functions
double A (double * y, double * z);
void DataBlocking(BaseFunction3D * p, MetroStep * trans_prob);

#endif
