#include <fstream>
#include <string>

using namespace std;

int main(){

	int n_props=4, n_nums;
	string prop[] = {"ene","heat","mag","chi"};
	double bin;

	ifstream input;
	ofstream output;

	for(int i=0; i<n_props; i++){
		n_nums=0;
		input.open("output/output."+prop[i]+".0");
		while(!input.eof()){
			input >> bin;
			n_nums++;
		}
		input.close();

		input.open("output/output."+prop[i]+".0");
		for(int num=0; num<n_nums-3; num++)
			input >> bin;
		output.open("data/data."+prop[i]+".dat",ios::app);
		input >> bin;
		output << bin << "	";
		input >> bin;
		output << bin << endl;
		output.close();
		input.close();
	}

	return 0;
}
		
		
