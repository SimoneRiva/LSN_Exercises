#ifndef _funzionebase_h_
#define _funzionebase_h_
#include <cmath>

class BaseFunction {
	public:
		virtual double Eval (double x) const =0;
		virtual void Set (double m, double s) =0;
};

class square_psi_trial : public BaseFunction {
	private:
		double mu, sigma;
	public:
		square_psi_trial(double m, double s) {mu=m; sigma=s;}
		virtual void Set(double m, double s) {mu=m; sigma=s;}
		double GetMu() {return mu;}
		double GetSigma() {return sigma;}
		virtual double Eval(double x) const {
			double psi = exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2.*sigma*sigma));
			return psi*psi;
		}
	};

#endif
