#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h"
#include "random.h"
#include "error.h"
#include "funzioneBase.h"

using namespace std;

int M=1e6, N=100, L=(M/N), k; //Total # of throws, # of blocks, throws per block
double sum;
double a=0., b=1., pmax=0.5*M_PI; //Integration interval and max of the distribution
vector<double> x;
vector<double> I (N,0.); //u=uniform sampling
vector<double> I2 (N,0.);
vector<double> sum_prog (N,0.);
vector<double> sum2_prog (N,0.);
vector<double> err_prog (N,0.);
vector<double> r;

void DataBlocking(BaseFunction * integ, BaseFunction * prob);

int main (){

//Setting random number generator
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
		   input >> property;
		   if( property == "RANDOMSEED" ){
		      input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		      rnd.SetRandom(seed,p1,p2);
		   }
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	BaseFunction * integ = new Integrand();
	BaseFunction * uniform = new UniformSampling();
	BaseFunction * importance = new ImportanceSampling();
	ofstream out;

	for (int h=0; h<N; h++)
		x.push_back(h);
//Sampling with a uniform distribution
	for (int h=0; h<M; h++)
		r.push_back(rnd.Rannyu(a, b));

	DataBlocking(integ, uniform);

	out.open("uniform_sampling.out");
	for (int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]-1.) << "	" << err_prog[h] << endl;
	out.close();

//sampling with p(x)=pi/2*(1-(1/2)*(pi*x/2)^2+(1/24)*(pi*x/2)^4)
	for (int h=0; h<M; h++)
		r[h]=(rnd.AccRej(importance, a, b, pmax));

	DataBlocking(integ, importance);

	out.open("importance_sampling.out");
	for (int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]-1.) << "	" << err_prog[h] << endl;
	out.close();

//This is to test the accept-reject method of Random: all the numbers produced are written in "accrejtest.out" and can be put in a histogrma to check if the expected distribution is obtained
/*
	out.open("accrejtest.out");
	for (int h=0; h<M; h++)
		out << r[h] << endl;
	out.close();
*/

	rnd.SaveSeed();
	return 0;
}


void DataBlocking(BaseFunction * integ, BaseFunction * prob){
	for(int i=0; i<N; i++){
		sum=0;
		for (int j=0; j<L; j++){
			k=j+i*L;
			double x = r[k];
			sum += (integ->Eval(x))/(prob->Eval(x));
		}
		I[i]=(b-a)*sum/L;
		I2[i]=I[i]*I[i];
		sum_prog[i]=0.;
		sum2_prog[i]=0.;
		for (int j=0; j<(i+1); j++){
			sum_prog[i] += I[j];
			sum2_prog[i] += I2[j];
		}
		sum_prog[i] /= (i+1);
		sum2_prog[i] /= (i+1);
		err_prog[i]=error(sum_prog, sum2_prog, i);
	}
}
