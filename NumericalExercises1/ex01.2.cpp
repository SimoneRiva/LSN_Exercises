#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

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

	int N[] = {1,2,10,100}; //Sizes of sets on which the average is computed
	double ave[3]; //Contains the averages for the three distributions
	double lambda=1., mu=0., gamma=1.;
	ofstream out1, out2, out3;

	out1.open("reg_dice.out");
	out2.open("exp_dice.out");
	out3.open("cauchy_dice.out");
	out1 << fixed << setprecision(6);
	out2 << fixed << setprecision(6);
	out3 << fixed << setprecision(6);

	for (int i=0; i<1e4; i++){ //10^4 realizations
		for (int k=0; k<4; k++){ //Index for N
			ave[0]=0;
			ave[1]=0;
			ave[2]=0;
			for (int j=0; j<N[k]; j++){
				ave[0] += rnd.Rannyu();
				ave[1] += rnd.RanExp(lambda);
				ave[2] += rnd.RanCauchy(mu, gamma);
			}
			ave[0]/=N[k];
			ave[1]/=N[k];
			ave[2]/=N[k];

			out1 << ave[0] << "	";
			out2 << ave[1] << "	";
			out3 << ave[2] << "	";
		}
		out1 << endl;
		out2 << endl;
		out3 << endl;
	}

	out1.close();
	out2.close();
	out3.close();

  rnd.SaveSeed();
	return 0;

}

/*

Generates 3 output files for regular, exponential and Cauchy dice, each organized as follows:

##N=1	N=2		N=10	N=100 -> this line is not present in the output files
S1		S2		S10		S100	-
S1		S2		S10		S100	 |
S1		S2		S10		S100	 |
...		...		...		...	   | 10^4 evaluations of the average SN of N numbers
...		...		...		...	   |
...		...		...		...	   |
S1		S2		S10		S100	-

*/

