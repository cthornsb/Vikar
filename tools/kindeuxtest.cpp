#include "vikar_core.h"

int main(){
	Kindeux kind;
	double exRecoil[2];
	exRecoil[0] = 0.7695;
	exRecoil[1] = 2.3200;
	kind.Initialize(7,2,8,1,-2.085,2,exRecoil,1.0);
	std::vector<std::string> fnames;
	fnames.push_back("/home/cory/Research/VANDLE/vikar2/xsections/ground_state.xsect");
	fnames.push_back("/home/cory/Research/VANDLE/vikar2/xsections/first_excited.xsect");
	fnames.push_back("/home/cory/Research/VANDLE/vikar2/xsections/second_excited.xsect");
	
	kind.SetDist(fnames, 16.0, 1E6);
	kind.Print();
	
	std::ofstream outFile("kind_test.dat");
	double my_array[3];
	for(unsigned int i = 0; i < 10000; i++){
		kind.Sample(my_array);
		outFile << my_array[0] << "\t" << my_array[1] << "\t" << my_array[2] << "\n";
	}
	outFile.close();
	
	return 0;
}
