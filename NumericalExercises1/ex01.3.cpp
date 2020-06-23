#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "random.h"
#include "error.h"

using namespace std;

double d=3., l=2.;
int M=1e6, N=100, Nthr=(M/N), k; //Total # of throws, # of blocks, throws per block
double Nhit; //Throws that intersect a line
vector<double> x; //x-coordinates
//Quantities for blocking average of pi
vector<double> pi (N,0.);
vector<double> pi2 (N,0.);
vector<double> sum_prog (N,0.);
vector<double> sum2_prog (N,0.);
vector<double> err_prog (N,0.);
vector<double> xend; //x-coordinates of one end of the needle
vector<double> theta; //angle formed with x-axis

int buffon (double x, double theta); // gives the results of the Buffon experiment with a needle of length l and a distance d between vertical lines and returns 1 if the needle intersects a line, 0 otherwise (gives a caution  message if d<l but still runs the experiment). The x-coordinate of one end of the needle and the angle measured form the horizontal direction must be given.
void DataBlockingPi();

int main (){

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


//Setting x-coordinates and generating random numbers
	for (int h=0; h<N; h++)
		x.push_back(h);
	for (int h=0; h<M; h++) {
		xend.push_back(rnd.Rannyu(0., d));
		theta.push_back(rnd.Rannyu(0., 2*M_PI));
	}

	DataBlockingPi();

	ofstream out;
	out.open("pi.out");
	for (int h=0; h<N; h++)
		out << x[h] << "	" << sum_prog[h] << "	" << err_prog[h] << endl;
	out.close();

   rnd.SaveSeed();
	return 0;
}


int buffon (double x, double theta){ // gives the results of the Buffon experiment with a needle of length l and a distance d between vertical lines and returns 1 if the needle intersects a line, 0 otherwise (gives a caution  message if d<l but still runs the experiment). The x-coordinate of one end of the needle and the angle measured form the horizontal direction must be given.
	if(d<=l) cout << "Caution: d is not greater than l !";

	double x2=x+l*cos(theta);
	if(floor(x/d)==floor(x2/d))
		return 0;
	else
		return 1;
}


void DataBlockingPi(){ //Performs the blocking averages of pi
	for(int i=0; i<N; i++){
		Nhit=0;
		for (int j=0; j<Nthr; j++){
			k=j+i*Nthr;
			Nhit += buffon(xend[k], theta[k]);
		}
		pi[i]=2.*l*Nthr/(Nhit*d);
		pi2[i]=pi[i]*pi[i];
		sum_prog[i]=0;
		sum2_prog[i]=0;
		for (int j=0; j<(i+1); j++){
			sum_prog[i] += pi[j];
			sum2_prog[i] += pi2[j];
		}
		sum_prog[i] /= (i+1);
		sum2_prog[i] /= (i+1);
		err_prog[i]=error(sum_prog, sum2_prog, i);
	}
}
