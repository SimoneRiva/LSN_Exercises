#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "error.h"
#include "ex05.h"

using namespace std;

int main (){

//Setting the random number generator
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

//|Psi|^2
  BaseFunction3D * s = new square_1s();
  BaseFunction3D * p = new square_2p();

//Transition probabilities
	MetroStep * uniform = new MetroStepUnif();
	MetroStep * gaussian = new MetroStepGauss();

// 1s //

	l=1.2;
	sigma=0.75;

//Uniform transition probability
  for(int i=0; i<3; i++)
    x[i]=10.;

	out_points.open("points/1s_points_unif.out");
	out_radii.open("radii/1s_radii_unif.out");

	DataBlocking(s, uniform);

	out_points.close();
	out_radii.close();

	cout << "The fraction of acceptance for the 1s orbital with uniform transition probability is " << double(accepted)/double(nstep) << endl;

//Gaussian transition probability
  for(int i=0; i<3; i++)
    x[i]=10.;

	out_points.open("points/1s_points_gauss.out");
	out_radii.open("radii/1s_radii_gauss.out");

	DataBlocking(s, gaussian);

	out_points.close();
	out_radii.close();

	cout << "The fraction of acceptance for the 1s orbital with gaussian transition probability is " << double(accepted)/double(nstep) << endl;

// 2p //

	l=3.;
	sigma=1.9;

//Uniform transition probability
  for(int i=0; i<3; i++)
    x[i]=10.;

	out_points.open("points/2p_points_unif.out");
	out_radii.open("radii/2p_radii_unif.out");

	DataBlocking(p, uniform);

	out_points.close();
	out_radii.close();

	cout << "The fraction of acceptance for the 2p orbital with uniform transition probability is " << double(accepted)/double(nstep) << endl;

//Gaussian transition probability
  for(int i=0; i<3; i++)
    x[i]=10.;

	out_points.open("points/2p_points_gauss.out");
	out_radii.open("radii/2p_radii_gauss.out");

	DataBlocking(p, gaussian);

	out_points.close();
	out_radii.close();

	cout << "The fraction of acceptance for the 2p orbital with gaussian transition probability is " << double(accepted)/double(nstep) << endl;

  rnd.SaveSeed();
	return 0;
}



double A (double * y, double * z, BaseFunction3D * p){
  return min(1., (p->Eval(y)/p->Eval(z)));
}

void DataBlocking(BaseFunction3D * p, MetroStep * trans_prob){
	step=0;
  while(step<equilstep){
    trans_prob->Step(p);
	}
	step=0;
	accepted=0;
  for(int iblock=0; iblock<nblocks; iblock++){
		sum=0.;
		ave_r[iblock]=0.;
		ave2_r[iblock]=0.;
		for(int i=0; i<blocksize; i++){
    	trans_prob->Step(p);
			sum += sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
			out_points << x[0] << "	" << x[1] << "	" << x[2] << endl;
		}
		r[iblock] = sum/blocksize;
		for(int jblock=0; jblock<=iblock; jblock++){
			ave_r[iblock] += r[jblock];
			ave2_r[iblock] += r[jblock]*r[jblock];
		}
		ave_r[iblock] /= (iblock+1);
		ave2_r[iblock] /= (iblock+1);
		err_r[iblock] = error(ave_r, ave2_r, iblock);

		out_radii << ave_r[iblock] << "	" << err_r[iblock] << endl;
	}
}

void MetroStepUnif::Step(BaseFunction3D * p) const {
  double y[3];
  for(int i=0; i<3; i++)
    y[i] = x[i]+l*rnd.Rannyu(-1,1);
  double r = rnd.Rannyu(), alpha = A(y, x, p);
  if(r<=alpha){
    for(int i=0;i<3; i++)
      x[i]=y[i];
		accepted++;
   }
  step++;
}

void MetroStepGauss::Step(BaseFunction3D * p) const {
  double y[3];
  for(int i=0; i<3; i++)
    y[i] = x[i]+rnd.Gauss(0., sigma);
  double r = rnd.Rannyu(), alpha = A(y, x, p);
  if(r<=alpha){
    for(int i=0;i<3; i++)
      x[i]=y[i];
		accepted++;
   }
  step++;
}
