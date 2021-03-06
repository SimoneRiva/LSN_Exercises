#ifndef _TSP_
#define _TSP_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

const int nstep = 20000; //Number of cooling steps
const int num = 200; //Number of steps at each beta
double beta [nstep];
int n[nstep];
int accepted, moves;

const int N = 32; //# of cities
double x[N], y[N]; //Cities positions
double length; //Path's length

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

Path path, newpath;

//Functions (see "TSP.cpp" for further explanations)
void Initialize();
void MetroStep(int step);
void SwapMutation();
void ShiftMutation();
void InversionMutation();

#endif
