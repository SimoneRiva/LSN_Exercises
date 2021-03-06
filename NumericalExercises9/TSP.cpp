#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

	int m1, m2, num;
	double ave_len;

	Initialize();

	ofstream output;
	output.open("length.out");

	for(int step=0; step<nstep; step++){

		FakeRoulette();

		m1 = Select();
		m2 = Select();

		parents[0]=paths[m1];
		parents[1]=paths[m2];
		paths[0]=parents[0];
		paths[1]=parents[1];

		Crossover();

		SwapMutation();
		ShiftMutation();
		InversionMutation();

		paths[0].check();
		paths[1].check();
		length[0]=paths[0].L1();
		length[1]=paths[1].L1();
		SortFirstTwo();

//Write best length and average length on the best half population on "length.out"
		ave_len = 0.;
		num = 0;
		for(int m=n/2; m<n; m++){
			ave_len += length[m];
			num++;
		}
		output << length[n-1] << "	" << ave_len/(double)num << endl;

		if((step+1)%1000==0)
			cout << "Step " << step+1 << " completed" << endl;
	}

	output.close();

	cout << endl << "Best path's lenght is " << paths[n-1].L1() << endl;

//Write the cities coordinates in the best order on "best_path.out"
	output.open("best_path.out");
	for(int i=0; i<N; i++)
		output << x[paths[n-1].get(i)-1] << "	" << y[paths[n-1].get(i)-1] << endl;
	output << x[paths[n-1].get(0)-1] << "	" << y[paths[n-1].get(0)-1] << endl;
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


//Set starting population
	for(int m=0; m<n; m++){
		Path Id;
		Id.check();
		paths[m].set(0,1);
		for(int i=1; i<N; i++){
			int c = (int)(rnd.Rannyu()*(Id.Np()-1)+1);
			paths[m].set(i,Id.get(c));
			Id.remove(c);
		}
		paths[m].check();
//Save each path's length
		length[m]=paths[m].L1();
	}

//Order starting population
	Sort();

}

//Swaps two paths
void SwapPaths(int m1, int m2){
	double a = length[m1];
	length[m1] = length[m2];
	length[m2] = a;
	Path * p = new Path();
	*p = paths[m1];
	paths[m1] = paths[m2];
	paths[m2] = *p;
}

//Orders paths with decreasing length
void Sort(){
	for(int m1=n; m1>1; m1--){
		for(int m2=0; m2<m1-1; m2++){
			if(length[m2]<length[m2+1])
				SwapPaths(m2, m2+1);
		}
	}
}

//Sorts the new generated paths to the right place
void SortFirstTwo(){
	int m;
	m=1;
	while(length[m]<length[m+1]){
		SwapPaths(m,m+1);
		m++;
	}
	m=0;
	while(length[m]<length[m+1]){
		SwapPaths(m,m+1);
		m++;
	}
}

//Sets the cumulative probabilities of being selected
void FakeRoulette(){
	for(int mi=0; mi<n; mi++){
		prob[mi]=0.;
		for(int mj=0; mj<mi+1; mj++)
			prob[mi]+=exp(-length[mj]);
	}
}

//Selects a parent with probability proportional to exp(-L1(x1,...,xN))
double Select(){
	double r = rnd.Rannyu(0., prob[n-1]);
	int m=0;
	while(r>prob[m])
		m++;
	return m;
}

//Generates a uniform random number i between 1 and N, cuts the paths at that position and completes them with the missing cities in the order in which they appear in the consort
void Crossover(){
	const double Pc = 0.6;

	double r = rnd.Rannyu();
	if(r<=Pc){
		int i = (int)rnd.Rannyu(1,N);
		int k;

		k=0;
		for(int j1=1; j1<N; j1++){
			int j2;
			for(j2=1; j2<i; j2++){
				if(parents[0].get(j1)==parents[1].get(j2)){
					j2=i+1;
				}
			}
			if(j2==i){
				paths[1].set(i+k, parents[0].get(j1));
				k++;
			}
		}

		k=0;
		for(int j1=1; j1<N; j1++){
			int j2;
			for(j2=1; j2<i; j2++){
				if(parents[1].get(j1)==parents[0].get(j2)){
					j2=i+1;
				}
			}
			if(j2==i){
				paths[0].set(i+k, parents[1].get(j1));
				k++;
			}
		}
	}

}

//For both sons chooses randomly two cities (except for the first) and swaps them
void SwapMutation(){
	const double Pm = 0.09;

	int i, j;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(1,N);
		paths[0].swap(i,j);
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(1,N);
		paths[1].swap(i,j);
	}
}

//For both sons chooses randomly a starting and an ending point for a block of cities and shifts them of a random (positive or negative) number of positions
void ShiftMutation(){
	const double Pm = 0.09;

	int i, j, k;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		k = (int)rnd.Rannyu(-i,N-j);
		for(int h=i; h<=j; h++){
			for(int hk=0; hk<abs(k); hk++)
				paths[0].swap(h,h+k/abs(k));
		}
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		k = (int)rnd.Rannyu(-i,N-j);
		for(int h=i; h<=j; h++){
			for(int hk=0; hk<abs(k); hk++)
				paths[1].swap(h,h+k/abs(k));
		}
	}
}

//For both sons chooses randomly a starting and an ending point for a block of cities and inverts them
void InversionMutation(){
	const double Pm = 0.09;

	int i, j;
	double r = rnd.Rannyu();
	if(r<=Pm){
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		for(int h=0; h<(j-i+1)/2; h++)
				paths[0].swap(i+h,j-h);
		i = (int)rnd.Rannyu(1,N);
		j = (int)rnd.Rannyu(i,N);
		for(int h=0; h<(j-i+1)/2; h++)
				paths[1].swap(i+h,j-h);
	}
}
