#include <iostream>
#include <fstream>

#include "vikar_core.h"
#include "materials.h"
#include "detectors.h"

double proper_value(const std::string &prompt_, const double &min_, bool ge_=false){
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

double proper_value(const std::string &prompt_, const double &min_, const double &max_, bool ge_=false){
	double output = -1;
	while(true){
		std::cout << " " << prompt_;
		std::cin >> output;
		if(!ge_){
			if(output > min_ && output < max_){ break; }
			std::cout << "  Error: Invalid value! Input must be > " << min_ << " and < " << max_ << ".\n";
		}
		else{
			if(output >= min_){ break; }
			std::cout << "  Error: Invalid value! Input must be >= " << min_ << " and <= " << max_ << ".\n";
		}
	}
	return output;
}

void help(){
	std::cout << "\n Type 'quit' at any time to exit...\n";
	std::cout << " Type 'help' to display this dialogue...\n\n";
	std::cout << " Select a value to input:\n";
	std::cout << "  (0) - Beta\n";
	std::cout << "  (1) - Gamma\n";
	std::cout << "  (2) - Total Energy (MeV)\n";
	std::cout << "  (3) - Total Energy (MeV/u)\n";
	std::cout << "  (4) - Kinetic Energy (MeV)\n";
	std::cout << "  (5) - Kinetic Energy (MeV/u)\n";
	std::cout << "  (6) - Momentum (MeV/c)\n";
	std::cout << "  (7) - Velocity (m/s)\n";
	std::cout << "  (8) - Brho (T*m)\n";
	std::cout << "  (9) - Erho (MJ/C)\n";
}

int main(int argc, char* argv[]){
	double Z, A, BE_A, v;
	double M_MeV, M_amu, M_kg;
	double beta, gamma, Brho;
	double total, kinetic;
	double momentum;
	
	Z = proper_value("Enter particle Z: ", 0.0);
	A = proper_value("Enter particle A: ", 0.0);
	BE_A = proper_value("Enter particle BE/A (MeV): ", 0.0, true);
	
	Particle part("", Z, A, BE_A);
	M_MeV = part.GetMass();
	M_amu = part.GetMassAMU();
	M_kg = part.GetMassKg();
	
	std::cout << "  Rest Mass: " << M_MeV << " MeV/c^2\n";
	std::cout << "  Rest Mass: " << M_amu << " amu\n";
	std::cout << "  Rest Mass: " << M_kg << " kg\n";	
	
	help();
	
	double value;
	std::string input;
	while(true){
		std::cout << "\n Enter selection (0-9): "; std::cin >> input;
		if(input == "quit"){ break; }
		else if(input == "help"){ 
			help();
			continue;
		}
		else if(input == "0"){ // Beta
			value = proper_value("Enter beta (0-1): ", 0.0, 1.0);
			v = value*c;
		}
		else if(input == "1"){ // Gamma
			value = proper_value("Enter gamma: ", 1.0);
			v = c*std::sqrt(1.0 - 1.0/(value*value));
		}
		else if(input == "2"){ // Total energy
			value = proper_value("Enter total energy (MeV): ", 0.0);
			v = part.GetVfromTE(value);
		}
		else if(input == "3"){ // Total energy per nucleon
			value = proper_value("Enter total energy (MeV/u): ", 0.0);
			v = part.GetVfromTE(value*M_amu);
		}
		else if(input == "4"){ // Kinetic energy
			value = proper_value("Enter kinetic energy (MeV): ", 0.0);
			v = part.GetVfromKE(value);
		}
		else if(input == "5"){ // Kinetic energy per nucleon
			value = proper_value("Enter kinetic energy (MeV/u): ", 0.0);
			v = part.GetVfromKE(value*M_amu);
		}
		else if(input == "6"){ // Momentum
			value = proper_value("Enter momentum (MeV/c): ", 0.0);
			v = part.GetVfromP(value);
		}
		else if(input == "7"){ // Velocity
			v = proper_value("Enter velocity (m/s): ", 0.0, c);
		}
		else if(input == "8"){ // Brho
			value = proper_value("Enter Brho (T*m): ", 0.0);
			v = part.GetVfromP(value*Z*e_charge*6.241506479963E12*c);
		}
		else if(input == "9"){ // Erho
			value = proper_value("Enter Erho (MJ/C): ", 0.0);
			double temp = std::pow((value*1E6)*6.241506479963E12*Z*e_charge*c*c/M_MeV, 2.0); // In units of c^4.
			v = std::sqrt(0.5*(std::sqrt(temp/(c*c) + 4*temp) - temp/(c*c)));
		}
		else{ 
			std::cout << "  Unrecognized input (" << input << ")\n";
			std::cout << "   Type 'help' for a list of options.\n";
			continue;
		}
		
		beta = part.GetBeta(v);
		gamma = part.GetGamma(v);
		total = part.GetTEfromV(v);
		kinetic = part.GetKEfromV(v);
		momentum = part.GetPfromV(v);
		Brho = gamma*part.GetMassKg()*v/part.GetCharge();

		std::cout << "  Beta: " << beta << std::endl;
		std::cout << "  Gamma: " << gamma << std::endl;
		std::cout << "  Total Energy: " << total << " MeV\n";
		std::cout << "  Total Energy: " << total/M_amu << " MeV/u\n";
		std::cout << "  Kinetic Energy: " << kinetic << " MeV\n";
		std::cout << "  Kinetic Energy: " << kinetic/M_amu << " MeV/u\n";
		std::cout << "  Momentum: " << momentum << " MeV/c\n";
		std::cout << "  Velocity: " << v << " m/s\n";
		std::cout << "  Brho: " << Brho << " Tm\n";
		std::cout << "  Erho: " << Brho*v/1E6 << " MJ/C\n";		
	}
	
	return 0;
}
