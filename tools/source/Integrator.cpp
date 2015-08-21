#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#define PI 3.1415926540

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [filename]\n";
}

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}
	
	std::ifstream inFile(argv[1]);
	if(!inFile.good()){
		std::cout << " Error: Failed to load input file\n";
		return 1;
	}
	
	// X should be angle in degrees
	// Y should be dsigma/domega in mb/sr
	std::vector<float> valx, valy;
	float fvalx, fvaly;
	unsigned int count = 0;
	
	// Load the data values
	while(true){
		inFile >> fvalx >> fvaly;
		if(inFile.eof()){ break; }
		
		valx.push_back(fvalx);
		valy.push_back(fvaly);
		count++;
	}
	inFile.close();

	std::cout << " Loaded " << count << " values from file\n";
	double integral = 0.0;
	
	// Integrate using the trapezoidal rule
	double x1, x2, y1, y2;
	for(unsigned int i = 0; i < valx.size()-1; i++){
		x1 = valx[i]*PI/180.0; y1 = valy[i]*std::sin(x1);
		x2 = valx[i+1]*PI/180.0; y2 = valy[i+1]*std::sin(x2);
		integral += 0.5*(x2-x1)*(y2+y1);
	}
	
	if(integral < 0.0){ integral *= -1; }
	std::cout << " The cross section is " << integral*2*3.14159 << " mb\n";
	
	return 0;
}
