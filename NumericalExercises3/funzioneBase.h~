#ifndef _funzionebase_h_
#define _funzionebase_h_
#include <cmath>

double S0=100; //Asset starting price
double T=1;	//Delivery time
double K=100; //Strike price
double r=0.1; //Risk-free interest rate
double sigma=0.25, sigma2=sigma*sigma; //Volatility

class BaseFunction {
	public:
		virtual double Eval (double x) const =0;
};

class CallPrice : public BaseFunction {
	public:
		virtual double Eval(double x) const {
			return exp(-r*T)*max(0., x-K);
		}
};

class PutPrice : public BaseFunction {
	public:
		virtual double Eval(double x) const {
			return exp(-r*T)*max(0., K-x);
		}
};

#endif
