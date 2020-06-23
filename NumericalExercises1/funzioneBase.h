#ifndef _funzionebase_h_
#define _funzionebase_h_

class BaseFunction {
	public:
		virtual double Eval (double x) const =0;
};

class id : public BaseFunction { //Identity function
	public:
		virtual double Eval(double x) const {
			return x;
		}
};

class squared_dev : public BaseFunction { //Returns the squared deviation from 1/2
	public:
		virtual double Eval(double x) const {
			return (x-0.5)*(x-0.5);
		}
};

#endif
