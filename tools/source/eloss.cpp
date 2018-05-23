#include <iostream>
#include <fstream>

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
	std::cout << "  SYNTAX: " << prog_name_ << " [matfile] {Z, A, q, BE/A, maxE, thick}\n";
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

	double Z, A, q, BE_A;
	double energy;
	double range;
	double eloss;
	double thickness;
	
	if(argc < 8){
		Z = proper_value("Enter particle Z: ");
		A = proper_value("Enter particle A: ");
		q = proper_value("Enter particle q: ");
		BE_A = proper_value("Enter particle BE/A (MeV): ", 0.0, true);
		energy = proper_value("Enter maximum energy (MeV): ", 0.0, false);
		thickness = proper_value("Enter material thickness (m): ", 0.0, false);
	}
	else{
		Z = strtod(argv[2], NULL);
		A = strtod(argv[3], NULL);
		q = strtod(argv[4], NULL);
		BE_A = strtod(argv[5], NULL);
		energy = strtod(argv[6], NULL);
		thickness = strtod(argv[7], NULL);
	}

	Particle part("part", Z, A, BE_A);

	double elow = part.GetKEfromV(0.02*3E8);

	RangeTable table;
	table.Init(1000, elow, energy, q, part.GetMass(), &mat);

	std::cout << "\n Type 'quit' to exit...\n";

	std::string input;
	while(true){
		std::cout << " Enter particle energy (MeV): "; std::cin >> input;
	
		if(input == "quit"){ break; }
	
		energy = strtod(input.c_str(), NULL);
		range = table.GetRange(energy);
		eloss = energy - table.GetEnergy(range - thickness);
		
		std::cout << "  Range = " << range << " m\n";
		std::cout << "  Eloss = " << eloss << " MeV\n";
		std::cout << "  Eremain = " << energy-eloss << " MeV\n";
	}
	
	return 0;
}
