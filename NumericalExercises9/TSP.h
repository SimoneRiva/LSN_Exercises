#ifndef _TSP_
#define _TSP_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

const int nstep = 20000; //Iterations

const int n = 100; //Population
const int N = 32; //# of cities
double x[N], y[N]; //Cities positions
double length [n]; //Paths lenght
double prob [n]; //Cumulative selection probabilty

//Permutations of [1,2,...,N]
class Path{
	private:
		int _Np; //# of cities
		int * _path;

	public:
		Path(); //Creates path [1,2,...,N]
		~Path();

		int Np(); //Returns Np
		int get(int i); //Returns i-th element
		void remove(int i); //Removes i-th element
		void set(int i, int a); //Sets i-th element to a
		void swap(int i, int j); //Swaps i-th and j-th elements
		void check(); //Checks if _path actually contains a permutation, if not gives an error message
		double L1(); //Evaluates the path's length

		Path& operator= (const Path& p);
};

Path * paths = new Path [n];
Path * parents = new Path [2];

//Functions (see "TSP.cpp" for further explanations)
void Initialize();
void SwapPaths(int m1, int m2);
void Sort();
void SortFirstTwo();
void FakeRoulette();
double Select();
void Crossover();
void SwapMutation();
void ShiftMutation();
void InversionMutation();

#endif
