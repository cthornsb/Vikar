#include <iostream>
#include <fstream>

#include "vikar_core.h"
#include "materials.h"
#include "detectors.h"

double proper_value(const std::string &prompt_, const double &min_=0.0, bool ge_=false){
	double output = -1;
	while(true){
		std::cout << " " << prompt_;
		std::cin >> output;
		if(!ge_){
			if(output > min_){ break; }
			std::cout << "  Error: Invalid value! Input must be > " << min_ << ".\n";
		}
		else{
			if(output >= min_){ break; }
			std::cout << "  Error: Invalid value! Input must be >= " << min_ << ".\n";
		}
	}
	return output;
}

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [matfile]\n";
}

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	Material mat;
	if(!mat.ReadMatFile(argv[1])){ 
		std::cout << " Failed to load material file '" << argv[1] << "'\n";
		return 1; 
	}

	double Z, A, BE_A;
	double start = 0.0;
	double thickness = 0.0;
	
	Z = proper_value("Enter particle Z: ");
	A = proper_value("Enter particle A: ");
	BE_A = proper_value("Enter particle BE/A (MeV): ", 0.0, true);
	start = proper_value("Enter particle start E (MeV): ", 0.0, false);
	thickness = proper_value("Enter material thickness (m): ", 0.0, false);

	Particle part("", Z, A, BE_A);
	part.SetMaterial(&mat, start, 0.1);
	
	double newE = part.GetTableNewE(start, thickness);

	std::cout << "  Eloss = " << start - newE << " MeV\n";
	std::cout << "  Efinal = " << newE << " MeV\n";
	
	return 0;
}
