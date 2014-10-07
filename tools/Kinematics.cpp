// angleConvert.cpp
// Cory Thornsberry
// April 25, 2014
// Convert an xmgrace x-section file from CoM angle to lab angle

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <time.h>

#include "vikar_core.h"

int main(int argc, char* argv[]){
	std::ifstream inFile;
	bool external_fname;
	if(argc >= 2){
		// Supply output filename
		external_fname = true;
	}
	else{
		// Use default filename
		external_fname = false;
	}
	
	Kindeux kind;	
	double Ebeam, Mbeam, Mtarg, Meject, Q;
	double start_angle, stop_angle;
	//double radius; // m
	//const double c = 3E8; // m/s
	//const double Mn = 939.565; // Mev/c^2
	unsigned short num_states;
	unsigned int num_trials, num_in_range;
	double *RecoilEx;
	std::cout << "  Enter beam Energy (MeV): "; std::cin >> Ebeam;
	std::cout << "  Enter beam Mass (A): "; std::cin >> Mbeam;
	std::cout << "  Enter target Mass (A): "; std::cin >> Mtarg;
	std::cout << "  Enter g.s. Q-value (MeV): "; std::cin >> Q;
	std::cout << "  Enter ejectile Mass (A): "; std::cin >> Meject;
	std::cout << "  Enter number of states: "; std::cin >> num_states;
	//std::cout << "  Enter radius of detectors: "; std::cin >> radius;
	if(num_states > 0){
		RecoilEx = new double[num_states];
		for(unsigned short i = 0; i < num_states; i++){
			std::cout << "   Enter recoil excitation for state " << i+1 << " (MeV): "; 
			std::cin >> RecoilEx[i];
		}
	}
	else{
		RecoilEx = new double[1];
		RecoilEx[0] = 0.0;
	}

	std::cout << "  Enter start angle (deg): "; std::cin >> start_angle; start_angle *= deg2rad;
	std::cout << "  Enter stop angle (deg): "; std::cin >> stop_angle; stop_angle *= deg2rad;
	std::cout << "  Enter number of trials: "; std::cin >> num_trials; num_in_range = 0;

	/*Ebeam = 30.3; Mbeam = 7; Mtarg = 2; Q = -2.085;
	Meject = 1; num_states = 3; RecoilEx = new double[num_states];
	RecoilEx[0] = 0.0; RecoilEx[1] = 0.7695; RecoilEx[2] = 2.3200;
	start_angle = 90.0; stop_angle = 180.0; num_trials = 100000; num_in_range = 0; radius = 0.5;*/

	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q, num_states, RecoilEx, 0.0);
	
	std::ofstream outFile;
	if(external_fname){ 
		outFile.open(argv[1]); 
		if(!outFile.good()){
			std::cout << "  Error: Failed to open the output file, check that the path is correct\n";
			delete[] RecoilEx;
			return 1;
		}
	}
	else{ outFile.open("Kinematics.out"); }
	
	double current_energy;
	//double ToF;
	Vector3 dummy_vector;
	while(num_in_range < num_trials){
		kind.FillVars(Ebeam, 0.0, 0.0, current_energy, dummy_vector);
		if(dummy_vector.axis[1] >= start_angle && dummy_vector.axis[1] <= stop_angle){ 
			outFile << to_str(dummy_vector.axis[1]*rad2deg) + "\t" + to_str(current_energy) << "\n"; 
			//ToF = radius*std::sqrt(kind.GetMeject()/(2*current_energy*1.60217657E-13*6.02214129E26));
			//outFile << to_str(dummy_vector.axis[1]*rad2deg) + "\t" + to_str(ToF) << "\n"; 
			num_in_range++;
		}
	}
	
	outFile.close();
	if(external_fname){ std::cout << " Wrote output file " << argv[1] << std::endl; }
	else{ std::cout << " Wrote output file Kinematics.out\n\n"; }
	
	delete[] RecoilEx;
	
	return 0;
}
