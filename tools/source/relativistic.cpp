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

int main(int argc, char* argv[]){
	double Z, A, BE_A, E, v;
	
	Z = proper_value("Enter particle Z: ");
	A = proper_value("Enter particle A: ");
	BE_A = proper_value("Enter particle BE/A (MeV): ", 0.0, true);
	E = proper_value("Enter kinetic energy (MeV): ");

	Particle part("", Z, A, BE_A);
	v = part.GetVfromKE(E);

	std::cout << "  Mass: " << part.GetMass() << " MeV/c^2\n";
	std::cout << "  Mass: " << part.GetMassAMU() << " amu\n";
	std::cout << "  Gamma: " << part.GetGamma(v) << std::endl;
	std::cout << "  Energy: " << part.GetKEfromV(v)/part.GetMassAMU() << " MeV/u\n";
	std::cout << "  Momentum: " << part.GetPfromV(v) << " MeV/c\n";
	std::cout << "  Velocity: " << v << " m/s\n";
	
	return 0;
}
