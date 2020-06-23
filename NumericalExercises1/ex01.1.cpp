#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "random.h"
#include "error.h"
#include "funzioneBase.h"

using namespace std;

int M=1e5; //Total throws
int N=100; //# of blocks
int L=(M/N); //Size of blocks
int k;
double sum;
vector<double> x;
//Blocking method quantities
vector<double> ave (N,0.);
vector<double> ave2 (N,0.);
vector<double> sum_prog (N,0.);
vector<double> sum2_prog (N,0.);
vector<double> err_prog (N,0.);
vector<double> r;

void DataBlocking(BaseFunction * func);

int main (){

//1.1.1

//Setting the random number generator
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

	ofstream out;

//Setting x-coordinates and generating random numbers
	for (int h=0; h<N; h++)
		x.push_back(h);
	for (int h=0; h<M; h++)
		r.push_back(rnd.Rannyu());

	BaseFunction * f1 = new id ();
	DataBlocking(f1);

	out.open("ave_r.out");
	for (int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]-0.5) << "	" << err_prog[h] << endl;
	out.close();

//1.1.2

	BaseFunction * f2 = new squared_dev ();
	DataBlocking(f2);

	out.open("var_r.out");
	for (int h=0; h<N; h++)
		out << x[h] << "	" << (sum_prog[h]-(1./12.)) << "	" << err_prog[h] << endl;
	out.close();

	x.clear();
	r.clear();

//1.1.3

	int intervals=100;
	double width = 1./intervals;
	int n=1e4; //# of throws for each j
	int jmax=100;
	double squarechi, rand_num;
	vector<int> num (intervals, 0); //# of throws in each interval

//Setting random numbers
	for (int h=0; h<n*jmax; h++)
		r.push_back(rnd.Rannyu());

	out.open("square_chi.out");

	for (int j=0; j<jmax; j++){
		squarechi=0.;
		for (int interv=0; interv<intervals; interv++)
			num[interv]=0;
		for (int k=0; k<n; k++){ //Searching in which interval each number falls
			rand_num = r[j*n+k];
			int interv = 0;
			while  (rand_num>(interv+1)*width)
				interv++;
			num[interv]++;
		}
		for (int interv=0; interv<intervals; interv++)
			squarechi += (num[interv]-n*width)*(num[interv]-n*width)/(n*width);
		out << squarechi << endl;
	}

	out.close();

	rnd.SaveSeed();
	return 0;
}


void DataBlocking(BaseFunction * func){ //Performs the blocking averages of the quantity func(r), where r are the random vector elements
	for(int i=0; i<N; i++){
		sum=0;
		for (int j=0; j<L; j++){
			k=j+i*L;
			sum += func->Eval(r[k]);
		}
		ave[i]=sum/L;
		ave2[i]=ave[i]*ave[i];
		sum_prog[i]=0;
		sum2_prog[i]=0;
		for (int j=0; j<(i+1); j++){
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] /= (i+1);
		sum2_prog[i] /= (i+1);
		err_prog[i]=error(sum_prog, sum2_prog, i);
	}
}
