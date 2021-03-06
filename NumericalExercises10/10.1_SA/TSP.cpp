#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

	Initialize();

	ofstream output;
	output.open("length.out");

	for(int step=0; step<nstep; step++){
		for(int m=0; m<n[step]; m++)
			MetroStep(step);
		output << beta[step] << "	" << length << endl;

		if((step+1)%100==0)
			cout << "Step " << step+1 << " completed" << endl;
	}

	output.close();

	cout << endl << "Best path's lenght is " << length << endl << "Acceptance rate: " << (double)accepted/(double)moves << endl;

	output.open("best_path.out");
	for(int i=0; i<N; i++)
		output << x[path.get(i)-1] << "	" << y[path.get(i)-1] << endl;
	output << x[path.get(0)-1] << "	" << y[path.get(0)-1] << endl;
	output.close();

	return 0;
}



Path::Path(){
	_Np=N;
	_path=new int [_Np];
	for(int i=0; i<_Np; i++) _path[i]=i+1;
}

Path::~Path(){}

int Path::Np(){
	return _Np;
}

int Path::get(int i){
	return _path[i];
}

void Path::remove(int i){
	for(int j=i; j<_Np-1; j++)	_path[j]=_path[j+1];
	_Np--;
}

void Path::set(int i, int a){
	_path[i]=a;
}

void Path::swap(int i, int j){
	int a=_path[i];
	_path[i]=_path[j];
	_path[j]=a;
}

void Path::check(){
	for(int i=0; i<_Np; i++){
		if(_path[i]<1) cout << "Error: " << i << "-th element is smaller than 1" << endl;
		if(_path[i]>_Np) cout << "Error: " << i << "-th element is bigger than " << this->Np() << endl;
		for(int j=0; j<i; j++){
			if(_path[i]==_path[j]) cout << "Error: " << i << "-th and " << j  << "-th elements are equal" << endl;
		}
	}
	if(_path[0]!=1) cout << "Error: 1st element is not 1" << endl;
}

Path& Path::operator=(const Path& p) {
	_Np=p._Np;
	delete [] _path;
	_path= new int [_Np];
	for (int i=0; i<_Np; i++)
		_path[i]=p._path[i];
	return *this;
}

double Path::L1(){
	double l = sqrt((x[_path[0]-1]-x[_path[N-1]-1])*(x[_path[0]-1]-x[_path[N-1]-1])+(y[_path[0]-1]-y[_path[N-1]-1])*(y[_path[0]-1]-y[_path[N-1]-1]));
	for (int i=0; i<N-1; i++)
		l += sqrt((x[_path[i]-1]-x[_path[i+1]-1])*(x[_path[i]-1]-x[_path[i+1]-1])+(y[_path[i]-1]-y[_path[i+1]-1])*(y[_path[i]-1]-y[_path[i+1]-1]));

	return l;
}


void Initialize(){
//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();

//Set the cities on a circumference
/*
	for(int i=0; i<N; i++){
		double theta=rnd.Rannyu(0., 2.*M_PI);
		x[i]=cos(theta);
		y[i]=sin(theta);
	}
*/

//Set the cities inside a square

	for(int i=0; i<N; i++){
		x[i]=rnd.Rannyu();
		y[i]=rnd.Rannyu();
	}


//Set starting path
	Path Id;
	Id.check();
	path.set(0,1);
	for(int i=1; i<N; i++){
		int c = (int)(rnd.Rannyu()*(Id.Np()-1)+1);
		path.set(i,Id.get(c));
		Id.remove(c);
	}
	path.check();
//Save path's length
	length=path.L1();

//Set annealing schedule
	for(int step=0; step<nstep; step++){
		beta[step] = (double)(step+1)/100.;
		n[step] = num;
	}

	accepted = 0;
}

void MetroStep(int step){
	newpath = path;

	SwapMutation();
	ShiftMutation();
	InversionMutation();

	newpath.check();
	double newlength = newpath.L1();
	double alpha = min(1., exp(-beta[step]*(newlength-length)));

	double r = rnd.Rannyu();
	if(r<alpha){
		path = newpath;
		length = newlength;
		accepted++;
	}
	moves++;
}

//Chooses randomly two cities (except for the first) and swaps them
void SwapMutation(){
	const double Pm = 0.4;

	int i, j;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(1,N);
		newpath.swap(i,j);
	}
}

//Chooses randomly a starting and an ending point for a block of cities and shifts them of a random (positive or negative) number of positions
void ShiftMutation(){
	const double Pm = 0.3;

	int i, j, k;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		k = (int)rnd.Rannyu(-i,N-j);
		for(int h=i; h<=j; h++){
			for(int hk=0; hk<abs(k); hk++)
				newpath.swap(h,h+k/abs(k));
		}
	}
}

//Chooses randomly a starting and an ending point for a block of cities and inverts them
void InversionMutation(){
	const double Pm = 0.3;

	int i, j;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		for(int h=0; h<(j-i+1)/2; h++)
				newpath.swap(i+h,j-h);
	}
}
