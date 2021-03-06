#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "MolDyn_NVE.h"
#include "error.h"

using namespace std;

int main(){

	Input(); //Inizialization
	int nconf = 1, k;

	for (int iblock=0; iblock<nblocks; ++iblock){
		for(int i_prop=0; i_prop<n_props; ++i_prop)
			sum_prop[i_prop]=0.;

		for(int istep=1; istep<=blocksize; ++istep){
			k=iblock*blocksize+istep;
			Move(); //Move particles with Verlet algorithm
			MeasureEachStep(); //Measures properties at each step and adds them to the variables sum_prop for block averaging
			if(k%iprint == 0) cout << "Number of time-steps: " << k << endl;
			if(k%10 == 0){
				Measure(); //Properties measurement
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}
			if(k == nstep-1) ConfOld(); //Write old configuration to restart
		}
		for(int bin=0; bin<nbins; ++bin){
			double r = bin_size*bin;
			double dv = 4./3.*M_PI*(pow(r+bin_size, 3)-pow(r,3));
			sum_prop[bin+igofr] /= (rho*(double)npart*dv);
		}

		for(int i_prop=0; i_prop<n_props; ++i_prop){
			ave[i_prop][iblock] = sum_prop[i_prop]/blocksize;
			sum_prog[i_prop][iblock]=0.;
			sum2_prog[i_prop][iblock]=0.;
			for(int jblock=0; jblock<(iblock+1); ++jblock){
				sum_prog[i_prop][iblock] += ave[i_prop][jblock];
				sum2_prog[i_prop][iblock] += ave[i_prop][jblock]*ave[i_prop][jblock];
			}
			sum_prog[i_prop][iblock] /= (iblock+1);
			sum2_prog[i_prop][iblock] /= (iblock+1);
			err[i_prop][iblock] = error(sum_prog[i_prop], sum2_prog[i_prop], iblock);
		}
	}

	ConfFinal(); //Write final configuration to restart

	ofstream out_etot, out_epot, out_ekin, out_temp, out_pres, out_gofr, out_gave;
	out_etot.open("averages/ave_etot.out");
	out_epot.open("averages/ave_epot.out");
	out_ekin.open("averages/ave_ekin.out");
	out_temp.open("averages/ave_temp.out");
	out_pres.open("averages/ave_pres.out");
	out_gofr.open("averages/gofr.out");
	out_gave.open("averages/gave.out");
	for(int iblock=0; iblock<nblocks; ++iblock){
		out_etot << sum_prog[ie][iblock] << "	" << err[ie][iblock] << endl;
		out_epot << sum_prog[iv][iblock] << "	" << err[iv][iblock] << endl;
		out_ekin << sum_prog[ik][iblock] << "	" << err[ik][iblock] << endl;
		out_temp << sum_prog[it][iblock] << "	" << err[it][iblock] << endl;
		out_pres << sum_prog[ip][iblock] << "	" << err[ip][iblock] << endl;
//g(r)
		for(int bin=0; bin<nbins; ++bin){
			out_gofr << iblock+1 << "	" << ave[bin+igofr][iblock] << endl;
			if(iblock==(nblocks-1))
				out_gave << sum_prog[bin+igofr][iblock] << "	" << err[bin+igofr][iblock] << endl;
		}
	}
	out_etot.close();
	out_epot.close();
	out_ekin.close();
	out_temp.close();
	out_pres.close();
	out_gofr.close();
	out_gave.close();

	return 0;
}


void Input(void){ //Prepare everything for the simulation
  ifstream ReadInput,ReadConf,ReadOldConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> old_config;
  if(old_config!=0 && old_config!=1) cerr << endl << endl << "Attention! The first line in 'input.dat' should contain 1 (if you want to read the configuration at time t-dt from 'old.0' or 0 (if you don't want to)" << endl << endl;
	ReadInput >> rescale_vel;

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

	blocksize = nstep/nblocks;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  ie = 0; //Total energy
  iv = 1; //Potential energy
  ik = 2; //Kinetic energy
  it = 3; //Temperature
	ip = 4; //Pressure
  n_props = 5; //Number of observables

//measurement of g(r)
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities if the old configuration is not provided
 if(old_config==0){
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
 }

//Read old configuration if provided
  if(old_config==1){
     double fs;
     cout << "Read old configuration from file old.0 " << endl << endl;
     ReadOldConf.open("old.0");
     for (int i=0; i<npart; ++i){
       ReadOldConf >> xold[i] >> yold[i] >> zold[i];
       xold[i] = xold[i] * box;
       yold[i] = yold[i] * box;
       zold[i] = zold[i] * box;
     }
     ReadOldConf.close();

//Correct r(t) to match the correct temperature if requested
     if(rescale_vel==1){
       Move();
//Temperature
       double t=0., stima_temp;
       for (int i=0; i<npart; ++i){
         vx[i] = (x[i]-xold[i])/delta;
         vy[i] = (y[i]-yold[i])/delta;
         vz[i] = (z[i]-zold[i])/delta;

         t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
       }
       stima_temp = (2.0/3.0) * t/(double)npart;
       fs = sqrt(temp/stima_temp);
//Rescale velocities and correct r(t)
       for (int i=0; i<npart; ++i){
         vx[i] *= fs;
         vy[i] *= fs;
         vz[i] *= fs;

         xold[i] = Pbc(x[i] - delta*vx[i]);
         yold[i] = Pbc(y[i] - delta*vy[i]);
         zold[i] = Pbc(z[i] - delta*vz[i]);
       }
     }
  }

  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("measures/output_epot.dat",ios::app);
  Ekin.open("measures/output_ekin.dat",ios::app);
  Temp.open("measures/output_temp.dat",ios::app);
  Etot.open("measures/output_etot.dat",ios::app);
  Pres.open("measures/output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = 48./pow(dr,12)-24./pow(dr,6);

//Potential energy
       v += vij;
       p += pij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_prop[iv] = v/(double)npart; //Potential energy per particle
    stima_prop[ik] = t/(double)npart; //Kinetic energy per particle
    stima_prop[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_prop[ie] = (t+v)/(double)npart; //Total energy per particle
    stima_prop[ip] = rho*stima_prop[it] + p/(3.0*vol); //Pressure

    Epot << stima_prop[iv]  << endl;
    Ekin << stima_prop[ik]  << endl;
    Temp << stima_prop[it] << endl;
    Etot << stima_prop[ie] << endl;
    Pres << stima_prop[ip] << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write old configuration
  ofstream WriteConf;

  cout << endl << "Print old configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){ //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}



void MeasureEachStep(){ //Measures properties at each step and adds them to the variables sum_prop for block averaging
  int bin;
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;
//reset the hystogram of g(r)
  for (int bin=igofr; bin<igofr+nbins; ++bin)
     stima_prop[bin]=0.;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
     for (int j=i+1; j<npart; ++j){

       dx = Pbc( x[i] - x[j] );
       dy = Pbc( y[i] - y[j] );
       dz = Pbc( z[i] - z[j] );

       dr = dx*dx + dy*dy + dz*dz;
       dr = sqrt(dr);

//update of the histogram of g(r)
			  if(dr < box/2.)
				  stima_prop[(int(dr/bin_size))+igofr]+=2.;

       if(dr < rcut){
         vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
         pij = 48./pow(dr,12)-24./pow(dr,6);

//Potential energy
         v += vij;
         p += pij;
       }
     }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i)
     t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  stima_prop[iv] = v/(double)npart; //Potential energy per particle
  stima_prop[ik] = t/(double)npart; //Kinetic energy per particle
  stima_prop[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_prop[ie] = (t+v)/(double)npart; //Total energy per particle
  stima_prop[ip] = rho*stima_prop[it] + p/(3*vol); //Pressure

	for(int i_prop=0; i_prop<n_props; ++i_prop)
     sum_prop[i_prop] += stima_prop[i_prop];

  return;
}
