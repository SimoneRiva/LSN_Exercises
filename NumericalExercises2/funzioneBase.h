#ifndef _funzionebase_h_
#define _funzionebase_h_

#include <iostream>
#include <cmath>

class BaseFunction {
public:
	virtual double Eval (double x) const =0;
};

class Integrand : public BaseFunction {
	public:
		virtual double Eval(double x) const {
			return 0.5*M_PI*cos(0.5*M_PI*x);
		}
};

class UniformSampling : public BaseFunction {
	public:
		virtual double Eval(double x) const {
			if(0<=x && x<=1)
				return 1.;
			else
				return 0.;
		}
};

class ImportanceSampling : public BaseFunction {
	public:
		virtual double Eval(double x) const {
			double norm=0.5*M_PI*(1.-(pow(M_PI, 2)/24.)+(pow(M_PI, 4)/1920.)); // normalization of p(x)
			return 0.5*M_PI*(1.-0.5*(0.5*M_PI*x)*(0.5*M_PI*x)+(pow((0.5*M_PI*x),4)/24.))/norm;
		}
};

#endif
