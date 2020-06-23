#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "VMC_1D.h"
#include "random.h"
#include "error.h"
#include "funzioneBase.h"

using namespace std;

int main (){

	Initialize();

	out_energy_search.open("energy_search.out");

	energy = TrialEnergy();
	for(int istep=0; istep<coolstep; istep++){
		acceptedSA=0;
		for(int m=0; m<n[istep]; m++)
			SAStep(istep);
//			cout << "The fraction of acceptance at step " << istep << " of the SA (beta=" << beta[istep] <<") is " << double(acceptedSA)/double(n[istep]) << endl;
			cout << "step " << istep+1 << " of " << coolstep << " completed" << endl;
	}

	out_energy_search.close();

	prob->Set(mu, sigma);
	BestEnergy();

  rnd.SaveSeed();

	return 0;
}


double A (double y, double z, BaseFunction * p){
  return min(1., (p->Eval(y)/p->Eval(z)));
}

void MetroStepUnif(BaseFunction * p){
  double y = x+l*rnd.Rannyu(-1,1);
  double r = rnd.Rannyu(), alpha = A(y, x, p);
  if(r<=alpha){
    x=y;
		accepted++;
   }
  step++;
}

void Initialize(){
//Set rnd
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

//Set annealing schedule
	for(int istep=0; istep<coolstep; istep++){
		beta[istep] = (double)(istep+1)*5.;
		n[istep] = num;
	}
}

double Eloc(double y){
	double e = -hbar*hbar/(2.*mass)*(((y-mu)*(y-mu)/pow(sigma,4)-1./(sigma*sigma))*exp(-(y-mu)*(y-mu)/(2.*sigma*sigma))+((y+mu)*(y+mu)/pow(sigma,4)-1./(sigma*sigma))*exp(-(y+mu)*(y+mu)/(2.*sigma*sigma)))/(exp(-(y-mu)*(y-mu)/(2.*sigma*sigma))+exp(-(y+mu)*(y+mu)/(2.*sigma*sigma)))+(pow(y,4)-2.5*y*y);
	return e;
}

//Performs a SA step by generating a couple of new parameters, evaluating the average energy and accepting the change with probability alpha (see below)
void SAStep(int istep){
	double old_mu = mu, old_sigma = sigma;
	mu += l_mu*rnd.Rannyu(-1, 1);
	sigma = max(0.1, sigma+l_sigma*rnd.Rannyu(-1, 1));

	prob->Set(mu, sigma);
	double new_energy = TrialEnergy();

	double alpha = min(1., exp(-beta[istep]*(new_energy-energy)));

	double r = rnd.Rannyu();
	if(r<alpha){
//		cout << energy << " -> " << new_energy << " accepted at beta = " << beta[istep] << endl;
		out_energy_search << fixed << setprecision(5) << beta[istep] << "	" << mu << "	" << sigma << "	" << new_energy << endl;
		energy = new_energy;
		acceptedSA++;
	}
	else{
		mu = old_mu;
		sigma = old_sigma;
	}
}

//Computes the block average energy with the current parameters and adds mu, sigma and ave_energy to "energy_search.out", returns the average for all blocks
double TrialEnergy(){ //Remember: always call prob->Set(mu, sigma); before calling TrialEnergy
	step=0;
	while(step<equilstep){ //This makes the Metropolis equilibrate
		MetroStepUnif(prob);
	}
	step=0;
	accepted=0;
	double ave_energy=0.;
  for(int iblock=0; iblock<nblocks; iblock++){
		sum=0.;
		for(int i=0; i<blocksize; i++){
   		MetroStepUnif(prob);
			sum += Eloc(x);
		}
		en[iblock] = sum/blocksize;
		ave_energy+=en[iblock]/nblocks;
	}
//	cout << "The fraction of acceptance in sampling psi square is " << double(accepted)/double(nstep) << endl;

	return ave_energy;
}

//Computes the block average energy with the current parameters mu and sigma and writes the points sampled by Metropolis in "points.out" and the block averages and errors in "energy.out"
void BestEnergy(){ //Remember: always call prob->Set(mu, sigma); before calling BestEnergy
	step=0;
	while(step<equilstep){ //This makes the Metropolis equilibrate
		MetroStepUnif(prob);
	}

	ofstream out_points, out_energy;
	out_points.open("points.out");
	out_energy.open("energy.out");
	step=0;
	accepted=0;
  for(int iblock=0; iblock<nblocks; iblock++){
		sum=0.;
		av_en[iblock]=0.;
		av2_en[iblock]=0.;
		for(int i=0; i<blocksize; i++){
   		MetroStepUnif(prob);
			sum += Eloc(x);
			out_points << x << endl;
		}
		en[iblock] = sum/blocksize;
		for(int jblock=0; jblock<=iblock; jblock++){
			av_en[iblock] += en[jblock];
			av2_en[iblock] += en[jblock]*en[jblock];
		}
		av_en[iblock] /= (iblock+1);
		av2_en[iblock] /= (iblock+1);
		err_en[iblock] = error(av_en, av2_en, iblock);

		out_energy << av_en[iblock] << "	" << err_en[iblock] << endl;
	}
	out_points.close();
	out_energy.close();

//	cout << "The fraction of acceptance in sampling psi square is " << double(accepted)/double(nstep) << endl;

	cout << endl << "The best values for the parameters are:" << endl << "mu = " << mu << endl << "sigma = " << sigma << endl << "The variational ground state energy is " << av_en[nblocks-1] << endl << endl;
}
