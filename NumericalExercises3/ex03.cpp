#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "random.h"
#include "error.h"
#include "funzioneBase.h"

using namespace std;

int M=1e6, N=100, L=(M/N), k; //Total # of throws, # of blocks, throws per block
double sum;
vector<double> x;
vector<double> Cave (N,0.);
vector<double> Cave2 (N,0.);
vector<double> sum_prog (N,0.);
vector<double> sum2_prog (N,0.);
vector<double> err_prog (N,0.);
vector<double> ST; //Asset final price

void DataBlocking(BaseFunction * func);

int main (){

// DIRECT SAMPLING

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

	BaseFunction * call = new CallPrice();
	BaseFunction * put = new PutPrice();
	ofstream out;

//Setting x-coordinates and generating the asset final price
	for(int h=0; h<N; h++)
		x.push_back(h);
	for(int h=0; h<M; h++)
		ST.push_back(S0*exp((r-0.5*sigma2)*T+sigma*rnd.Gauss(0., T)));

	DataBlocking(call);

	out.open("call_direct_sampling.out");
	for(int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]) << "	" << err_prog[h] << endl;
	out.close();

	DataBlocking(put);

	out.open("put_direct_sampling.out");
	for(int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]) << "	" << err_prog[h] << endl;
	out.close();


// DISCRETE STEPS SAMPLING

	int Nsteps=100;
	double Sti; //S at each step
	double tau=T/(double)Nsteps; //Step width

//Generating the GBM for the asset price step by step
	for(int h=0; h<M; h++){
		Sti=S0;
		for(int step=0; step<Nsteps; step++)
			Sti=Sti*exp((r-0.5*sigma2)*tau+sigma*rnd.Gauss(0., 1.)*sqrt(tau));
		ST[h]=Sti;
	}

	DataBlocking(call);

	out.open("call_discrete_sampling.out");
	for(int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]) << "	" << err_prog[h] << endl;
	out.close();

	DataBlocking(put);

	out.open("put_discrete_sampling.out");
	for(int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]) << "	" << err_prog[h] << endl;
	out.close();

	rnd.SaveSeed();
	return 0;
}


void DataBlocking(BaseFunction * func){ //Performs the blocking averages of the quantity func(r), where r are the random vector elements
	for(int i=0; i<N; i++){
		sum=0;
		for (int j=0; j<L; j++){
			k=j+i*L;
			sum += func->Eval(ST[k]);
		}
		Cave[i]=sum/L;
		Cave2[i]=Cave[i]*Cave[i];
		sum_prog[i]=0;
		sum2_prog[i]=0;
		for (int j=0; j<(i+1); j++){
			sum_prog[i] += Cave[j];
			sum2_prog[i] += Cave2[j];
		}
		sum_prog[i] /= (i+1);
		sum2_prog[i] /= (i+1);
		err_prog[i]=error(sum_prog, sum2_prog, i);
	}
}
