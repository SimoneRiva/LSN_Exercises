#include <cmath>
#include <vector>
#include "error.h"

using namespace std;

double error(vector<double> av, vector<double> av2, int n){

	if (n==0)
		return 0;
	else
		return sqrt((av2[n] - av[n]*av[n])/n);
}
