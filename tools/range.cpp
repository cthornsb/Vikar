#include <iostream>
#include <fstream>

#include "vikar_core.h"
#include "materials.h"
#include "detectors.h"
#include "structures.h"

int main(int argc, char *argv[]){
	if(argc < 7){
		std::cout << " Missing a required argument. Expected 6, received " << argc-1 << std::endl;
		std::cout << "  SYNTAX: ./range {material_filename} {particle_Z} {particle_A} {start_E} {stop_E} {points} [BE/A]\n";
		return 1;
	}
	
	double Z = (double)atof(argv[2]);
	double A = (double)atof(argv[3]);
	double BE_A = 0.0;

	double start = (double)atof(argv[4]);
	double stop = (double)atof(argv[5]);
	unsigned int points = (unsigned int)atoi(argv[6]);
	double step = (stop-start)/points;

	if(argc >= 8){
		BE_A = (double)atof(argv[7]);
	}
	double mass = (Z*proton_RME+(A-Z)*neutron_RME)-BE_A*A;
	
	std::cout << " Using particle with Z = " << Z << " and A = " << A << std::endl;
	std::cout << " Rest Mass: " << mass << " MeV\n";
	std::cout << " Start E: " << start << " MeV\n";
	std::cout << " Stop E: " << stop << " MeV\n";
	std::cout << " Step: " << step << " MeV\n";	
	
	Material mat;
	if(!mat.ReadMatFile(argv[1])){ 
		std::cout << " Failed to load material file '" << argv[1] << "'\n";
		return 1; 
	}
	
	std::cout << " Successfully loaded material '" << mat.GetName() << "'\n";
	mat.Print();

	std::ofstream output("range.out");
	for(unsigned int i = 0; i < points; i++){
		output << step*(i+1) << "\t" << mat.Range(step*(i+1), Z, mass) << "\n";
	}
	output.close();
	std::cout << " Done\n";
	
	return 0;
}
