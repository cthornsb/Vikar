// Kinematics.cpp
// Cory Thornsberry
// April 25, 2014
// Calculate the kinematics for an input reaction

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <time.h>

#include "vikar_core.h"
#include "kindeux.h"

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
	unsigned short num_states;
	unsigned int num_trials, num_in_range;
	double *RecoilEx;
	std::cout << "  Enter beam Energy (MeV): "; std::cin >> Ebeam;
	std::cout << "  Enter beam Mass (A): "; std::cin >> Mbeam;
	std::cout << "  Enter target Mass (A): "; std::cin >> Mtarg;
	std::cout << "  Enter g.s. Q-value (MeV): "; std::cin >> Q;
	std::cout << "  Enter ejectile Mass (A): "; std::cin >> Meject;
	std::cout << "  Enter number of states: "; std::cin >> num_states;

	// Declare storage arrays
	if(num_states == 0){ num_states = 1; }
	RecoilEx = new double[num_states];
	for(unsigned int i = 0; i < num_states; i++){
		if(i > 0){
			std::cout << "   Enter recoil excitation for state " << i << " (MeV): "; 
			std::cin >> RecoilEx[i];
		}
		else{
			std::cout << "   Using recoil excitation of 0.0 MeV for g.s.\n";
			RecoilEx[0] = 0.0;
		}
	}

	std::cout << "  Enter start lab angle (deg): "; std::cin >> start_angle; start_angle *= deg2rad;
	std::cout << "  Enter stop lab angle (deg): "; std::cin >> stop_angle; stop_angle *= deg2rad;
	std::cout << "  Enter number of trials: "; std::cin >> num_trials; num_in_range = 0;

	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q, num_states, RecoilEx);
	
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
	
	outFile << "Kinematics\toutput file v. 1.0\n";
	outFile << "Ebeam\t" << Ebeam << " MeV\n";
	outFile << "Mbeam\t" << Mbeam << " AMU\n";
	outFile << "Mtarg\t" << Mtarg << " AMU\n";
	outFile << "Meject\t" << Meject << " AMU\n";
	outFile << "num_states\t" << num_states << "\n";
	for(unsigned int i = 0; i < num_states; i++){
		outFile << "RecoilEx[" << i << "]\t" << RecoilEx[i] << " MeV\n";
	}
	outFile << "start_angle\t" << start_angle*rad2deg << " degrees (LAB frame)\n";
	outFile << "stop_angle\t" << stop_angle*rad2deg << " degrees (LAB frame)\n";
	outFile << "num_trials\t" << num_trials << "\n";
	outFile << "CoMAngle\tEjectAngle\tEjectE\tRecoilAngle\tRecoilE\n";
	
	Vector3 eject, recoil;
	reactData rdata;
	while(num_in_range < num_trials){
		kind.FillVars(rdata, eject, recoil);
		if(eject.axis[1] >= start_angle && eject.axis[1] <= stop_angle){ 
			outFile << rdata.comAngle*rad2deg << "\t" << eject.axis[1]*rad2deg << "\t" << rdata.Eeject << "\t" << recoil.axis[1]*rad2deg << "\t" << rdata.Erecoil << "\n"; 
			num_in_range++;
		}
	}
	
	outFile.close();
	if(external_fname){ std::cout << " Wrote output file " << argv[1] << std::endl; }
	else{ std::cout << " Wrote output file Kinematics.out\n\n"; }
	
	delete[] RecoilEx;
	
	return 0;
}
