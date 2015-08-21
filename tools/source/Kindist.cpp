// Kindist.cpp
// Cory Thornsberry
// Oct. 2nd, 2014
// Calculate kinematics curves for given reactions

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <time.h>

#include "vikar_core.h"

int main(int argc, char* argv[]){
	Kindeux kind;	
	double Ebeam, Mbeam, Mtarg, Meject, Q;
	double start_angle, stop_angle;
	unsigned int num_states;
	unsigned int num_steps;
	double *RecoilEx;
	double Eeject, Erecoil;
	double Eeject2, Erecoil2;
	std::cout << "  Enter beam Energy (MeV): "; std::cin >> Ebeam;
	std::cout << "  Enter beam Mass (A): "; std::cin >> Mbeam;
	std::cout << "  Enter target Mass (A): "; std::cin >> Mtarg;
	std::cout << "  Enter g.s. Q-value (MeV): "; std::cin >> Q;
	std::cout << "  Enter ejectile Mass (A): "; std::cin >> Meject;
	std::cout << "  Enter number of states: "; std::cin >> num_states;

	// Open the output file
	std::ofstream outFile;
	if(argc >= 2){ 
		outFile.open(argv[1]); 
		if(!outFile.good()){
			std::cout << "  Error: Failed to open the output file, check that the path is correct\n";
			return 1;
		}
	}
	else{ outFile.open("Kindist.out"); }

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

	std::cout << "  Enter start com angle (deg): "; std::cin >> start_angle; start_angle *= deg2rad;
	std::cout << "  Enter stop com angle (deg): "; std::cin >> stop_angle; stop_angle *= deg2rad;
	std::cout << "  Enter number of steps: "; std::cin >> num_steps;

	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q, num_states, RecoilEx);
	
	outFile << "Kindist\toutput file v. 1.0\n";
	outFile << "Ebeam\t" << Ebeam << " MeV\n";
	outFile << "Mbeam\t" << Mbeam << " AMU\n";
	outFile << "Mtarg\t" << Mtarg << " AMU\n";
	outFile << "Meject\t" << Meject << " AMU\n";
	outFile << "num_states\t" << num_states << "\n";
	for(unsigned int i = 0; i < num_states; i++){
		outFile << "RecoilEx[" << i << "]\t" << RecoilEx[i] << " MeV\n";
	}
	outFile << "start_angle\t" << start_angle*rad2deg << " degrees (COM frame)\n";
	outFile << "stop_angle\t" << stop_angle*rad2deg << " degrees (COM frame)\n";
	outFile << "num_steps\t" << num_steps << "\n";
	outFile << "CoMAngle\tEjectTheta1\tEjectE1\tEjectTheta2\tEjectE2\tRecoilTheta1\tRecoilE1\tRecoilTheta2\tRecoilE2\n";
	
	Vector3 eject, recoil, eject2, recoil2;	
	double com_angle;
	double step_size = (stop_angle - start_angle)/num_steps;
	for(unsigned int i = 0; i < num_steps; i++){
		outFile << (start_angle+i*step_size)*rad2deg;
		for(int j = 0; j < (int)num_states; j++){
			kind.FillVars(Ebeam, Eeject, Erecoil, eject, recoil, com_angle, j, 0, (start_angle+i*step_size));
			kind.FillVars(Ebeam, Eeject2, Erecoil2, eject2, recoil2, com_angle, j, 1, (start_angle+i*step_size));
			
			outFile << "\t" << eject.axis[1]*rad2deg << "\t" << Eeject << "\t" << eject2.axis[1]*rad2deg << "\t" << Eeject2 << "\t"; // Ejectile
			outFile << recoil.axis[1]*rad2deg << "\t" << Erecoil << "\t" << recoil2.axis[1]*rad2deg << "\t" << Erecoil2; // Recoil
		}
		outFile << std::endl;
	}
	
	outFile.close();
	if(argc >= 2){ std::cout << " Wrote output file " << argv[1] << std::endl; }
	else{ std::cout << " Wrote output file Kindist.out\n\n"; }
	
	delete[] RecoilEx;
	
	return 0;
}
