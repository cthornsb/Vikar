#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

int main(int argc, char* argv[]){
	if(argc < 2){ 
		std::cout << " Error: No filename specified\n";
		return 1; 
	}
	
	std::ifstream inFile(argv[1]);
	if(!inFile.good()){
		std::cout << " Error: Failed to load input file\n";
		return 1;
	}
	
	std::vector<float> valx, valy;
	float fvalx, fvaly;
	
	// Load the data values
	while(true){
		inFile >> fvalx >> fvaly;
		if(inFile.eof()){ break; }
		
		valx.push_back(fvalx);
		valy.push_back(fvaly);
	}
	inFile.close();
	
	unsigned int count = valx.size();
	if(valx.size() < valy.size()){
		std::cout << " Warning! Found more Y values than X values\n";
		count = valx.size();
	}
	else if(valx.size() > valy.size()){
		std::cout << " Warning! Found more X values than Y values\n";
		count = valy.size();
	}
	
	std::cout << " Loaded " << count << " values from file\n";
	double integral = 0.0;
	
	// Integrate using the trapezoidal rule
	for(unsigned int i = 0; i < count-1; i++){
		integral += (valx[i+1]-valx[i])*(valy[i]*std::sin(valx[i])+valy[i+1]*std::sin(valx[i+1]))/2.0;
	}
	
	std::cout << " The integral is " << integral*2*3.14156 << std::endl;
	
	return 0;
}
