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
	std::cout << "  SYNTAX: " << prog_name_ << " [filename]\n";
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
	double start, stop;
	unsigned int points;
	
	Z = proper_value("Enter particle Z: ");
	A = proper_value("Enter particle A: ");
	stop = proper_value("Enter particle final E: ", 0.0, true);
	start = proper_value("Enter particle start E: ", stop, false);
	BE_A = proper_value("Enter particle BE/A: ", 0.0, true);
	std::cout << " Enter number points: "; std::cin >> points;
	
	double step = (start-stop)/points;
	double mass = (Z*proton_RME+(A-Z)*neutron_RME)-BE_A*A;

	std::ofstream output("range.out");
	
	output << "Z\t" << Z << "\n";
	output << "A\t" << A << "\n";
	output << "RestMass\t" << mass << " MeV\n";
	output << "StartE\t" << start << " MeV\n";
	output << "StopE\t" << stop << " MeV\n";
	output << "Step\t" << step << " MeV\n";	
	output << "Material\t" << mat.GetName() << "\n";
	output << "Filename\t" << argv[1] << "\n";
	mat.Print(&output);
	output << "Energy\tRange\n";

	std::cout << " Processing " << points << " data points... \n";
	for(unsigned int i = 0; i < points; i++){
		output << step*(i+1) << "\t" << mat.Range(step*(i+1), Z, mass) << "\n";
	}
	output.close();
	std::cout << "done\n";
	
	return 0;
}
