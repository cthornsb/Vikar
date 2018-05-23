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
	std::cout << "  SYNTAX: " << prog_name_ << " [matfile] <options>\n";
	std::cout << "   Available options:\n";
	std::cout << "    --table | Output range table 'range.out'.\n";
}

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	bool make_table = false;
	int index = 2;
	while(index < argc){
		if(strcmp(argv[index], "--table") == 0){
			std::cout << " Saving range table to 'range.out'.\n";
			make_table = true;
		}
		else{ 
			std::cout << " Error! Unrecognized option '" << argv[index] << "'!\n";
			help(argv[0]);
			return 1;
		}
		index++;
	}

	Material mat;
	if(!mat.ReadMatFile(argv[1])){ 
		std::cout << " Failed to load material file '" << argv[1] << "'\n";
		return 1; 
	}
	mat.Print();

	double Z, A, BE_A;
	double start = 0.0, stop;
	unsigned int points;
	
	Z = proper_value("Enter particle Z: ");
	A = proper_value("Enter particle A: ");
	BE_A = proper_value("Enter particle BE/A: ", 0.0, true);
	if(!make_table){ stop = proper_value("Enter particle E: ", 0.0, true); }
	else{
		stop = proper_value("Enter particle final E: ", 0.0, true);
		start = proper_value("Enter particle start E: ", stop, false);
		std::cout << " Enter number points: "; std::cin >> points;
	}
	
	double step = (start-stop)/points;
	double mass = (Z*proton_RME+(A-Z)*neutron_RME)-BE_A*A;

	if(!make_table){
		std::cout << "  Range = " << mat.Range(stop, Z, mass) << " m\n";
	}
	else{
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
	}
	
	return 0;
}
