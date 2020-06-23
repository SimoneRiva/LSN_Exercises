#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "random.h"
#include "ranWalk.h"

using namespace std;

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


	int M=1e4, N=100, L=(M/N), k; //Total # of simulations, # of blocks, simulations per block
	int len=100; //Random walk length (i.e. number of steps)
	double a=1.; //Random walk step
	double l_d, l_c, ave_d, ave2_d, ave_c, ave2_c, err_d, err_c; //Quantities for blocking averages
	ofstream out_d, out_c;

	Vector3D ** RW_d = new Vector3D * [M]; //d=discrete random walk
	Vector3D ** RW_c = new Vector3D * [M]; //c=continuum random walk
	for(int h=0; h<M; h++){
		RW_d[h] = new Vector3D ();
		RW_c[h] = new Vector3D ();
	}

	out_d.open("discrete_random_walk.out");
	out_c.open("continuum_random_walk.out");

	for(int step=1; step<len; step++){ //At each step
		for(int h=0; h<M; h++){ //Update all simulations of both RWs
			RW_d[h]->discRWStep(a, rnd.Rannyu(0., 3.), rnd.Rannyu(0., 2.));
			RW_c[h]->contRWStep(a, acos(1.-2.*rnd.Rannyu()), rnd.Rannyu(0., 2.*M_PI));
		}
//Compute the averages
		ave_d=0.;
		ave2_d=0.;
		ave_c=0.;
		ave2_c=0.;
		for(int i=0; i<N; i++){
			l_d=0;
			l_c=0;
			for (int j=0; j<L; j++){
				k=j+i*L;
				l_d += RW_d[k]->squareR();
				l_c += RW_c[k]->squareR();
			}
			ave2_d+=l_d/L;
			ave_d+=sqrt(l_d/L);
			ave2_c+=l_c/L;
			ave_c+=sqrt(l_c/L);
		}

		ave_d /= N;
		ave2_d /= N;
		err_d=sqrt((ave2_d-ave_d*ave_d)/N);
		ave_c /= N;
		ave2_c /= N;
		err_c=sqrt((ave2_c-ave_c*ave_c)/N);

		out_d << step << "	" << ave_d << "	" << err_d << endl;
		out_c << step << "	" << ave_c << "	" << err_c << endl;
	}

	out_d.close();
	out_c.close();

	rnd.SaveSeed();
	return 0;
}
