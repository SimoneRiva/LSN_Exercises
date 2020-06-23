#ifndef _VMC_1D_
#define _VMC_1D_

#include <fstream>
#include "random.h"
#include "funzioneBase.h"

using namespace std;

Random rnd;

const int coolstep = 100; //Number of cooling steps
const int num = 100; //Number of MC steps at each beta
double beta [coolstep];
int n[coolstep];
int acceptedSA=0; //New values for the parameters accepted by the SA at each beta

const int nstep=1e5, equilstep=1e4, nblocks=100, blocksize=nstep/nblocks; //Data blocking constants
int step=0, accepted=0; //Metropolis total and accepted moves at each trial

double l=1.2; //Max movement for the Metropolis
double mu=1., sigma=0.5, energy;
double l_mu=0.5, l_sigma=0.2; //Max parameters variations for the SA
double x=1.; //Starting point
double hbar=1., mass=1.;

BaseFunction * prob = new square_psi_trial(mu, sigma);
double en [nblocks], av_en [nblocks], av2_en [nblocks], err_en [nblocks];
double sum;

ofstream out_energy_search;

//Functions
double T (double y, double z);
double A (double y, double z);
void MetroStepUnif(BaseFunction * p);

void Initialize();
double Eloc(double y);
void SAStep(int istep);
double TrialEnergy();
void BestEnergy();

#endif
