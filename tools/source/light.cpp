#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "vandmc_core.hpp"
#include "materials.hpp"
#include "detectors.hpp"

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
	std::cout << "  SYNTAX: " << prog_name_ << " [matfile] {Z, A, BE/A, maxE, L0, kB, C}\n";
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
	
	mat.Print();	

	double Z, A, BE_A;
	double energy;

	// Birks' constants.
	double L0, kB, C;
	
	if(argc < 9){
		Z = proper_value("Enter particle Z: ");
		A = proper_value("Enter particle A: ");
		BE_A = proper_value("Enter particle BE/A (MeV): ", 0.0, true);
		energy = proper_value("Enter maximum energy (MeV): ", 0.0, false);
		L0 = proper_value("Enter L0 for Birks' (1/MeV): ", 0.0, false);
		std::cout << " Enter kB for Birks' (m/MeV): "; std::cin >> kB;
		std::cout << " Enter C for Birks' (m/MeV)^2: "; std::cin >> C;
	}
	else{
		Z = strtod(argv[2], NULL);
		A = strtod(argv[3], NULL);
		BE_A = strtod(argv[4], NULL);
		energy = strtod(argv[5], NULL);
		L0 = strtod(argv[6], NULL);
		kB = strtod(argv[7], NULL);
		C = strtod(argv[8], NULL);
	}

	Particle part("part", Z, A, BE_A);

	double elow = part.GetKEfromV(0.02*3E8);

	RangeTable table;
	table.Init(100, elow, energy, Z, part.GetMass(), &mat);
	table.InitBirks(L0, kB, C);
	
	table.Print();

	std::cout << "\n Type 'quit' to exit...\n";

	std::string input;
	while(true){
		std::cout << " Enter particle energy (MeV/ADC): "; std::cin >> input;
	
		if(input == "quit"){ break; }
	
		energy = strtod(input.c_str(), NULL);
		
		//std::cout << "  Range = " << table.GetRange(energy) << " m\n";
		std::cout << "  Light = " << table.GetKEfromLR(energy) << " MeV\n";
	}
	
	return 0;
}
