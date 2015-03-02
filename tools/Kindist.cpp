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

	std::cout << "  Enter start CoM angle (deg): "; std::cin >> start_angle; start_angle *= deg2rad;
	std::cout << "  Enter stop CoM angle (deg): "; std::cin >> stop_angle; stop_angle *= deg2rad;
	std::cout << "  Enter number of steps: "; std::cin >> num_steps;

	/*Ebeam = 30.3; Mbeam = 7; Mtarg = 2; Q = -2.085;
	Meject = 1; num_states = 3; RecoilEx = new double[num_states];
	RecoilEx[0] = 0.0; RecoilEx[1] = 0.7695; RecoilEx[2] = 2.3200;
	EjectTheta = new double[num_states];
	RecoilTheta = new double[num_states];
	Eeject = new double[num_states];
	Erecoil = new double[num_states];
	RecoilTheta2 = new double[num_states];
	Eeject2 = new double[num_states];
	Erecoil2 = new double[num_states];
	for(unsigned int i = 0; i < num_states; i++){ 
		EjectTheta[i] = 0.0; RecoilTheta[i] = 0.0; Eeject[i] = 0.0; Erecoil[i] = 0.0; 
		RecoilTheta[i] = 0.0; Eeject[i] = 0.0; Erecoil[i] = 0.0;
	}
	start_angle = 0.0; stop_angle = pi; num_steps = 180;*/

	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q, num_states, RecoilEx, 0.0);
	
	outFile << "# Kindist kinematics output file v. 1.0\n";
	outFile << "# Ebeam = " << Ebeam << " MeV\n";
	outFile << "# Mbeam = " << Mbeam << " AMU\n";
	outFile << "# Mtarg = " << Mtarg << " AMU\n";
	outFile << "# Meject = " << Meject << " AMU\n";
	outFile << "# num_states = " << num_states << "\n";
	for(unsigned int i = 0; i < num_states; i++){
		outFile << "#  RecoilEx[" << i << "] = " << RecoilEx[i] << " MeV\n";
	}
	outFile << "# start_angle = " << start_angle*rad2deg << " degrees (CoM frame)\n";
	outFile << "# stop_angle = " << stop_angle*rad2deg << " degrees (CoM frame)\n";
	outFile << "# num_steps = " << num_steps << "\n";
	outFile << "#################################################################################################################################\n";
	outFile << "# CoMAngle----EjectThetaA----EjectE1A----EjectE2A----RecoilTheta1A----RecoilE1A----RecoilTheta2A----RecoilE2A----EjectThetaB... #\n";
	
	Vector3 eject, recoil, eject2, recoil2;	
	double step_size = (stop_angle - start_angle)/num_steps;
	for(unsigned int i = 0; i < num_steps; i++){
		outFile << (start_angle+i*step_size);
		for(int j = 0; j < (int)num_states; j++){
			kind.FillVars(Ebeam, Eeject, Erecoil, eject, recoil, j, 0, (start_angle+i*step_size));
			kind.FillVars(Ebeam, Eeject2, Erecoil2, eject2, recoil2, j, 1, (start_angle+i*step_size));
			
			outFile << "\t" << eject.axis[1] << "\t" << Eeject << "\t" << Eeject2 << "\t"; // Ejectile
			outFile << recoil.axis[1] << "\t" << Erecoil << "\t" << recoil2.axis[1] << "\t" << Erecoil2; // Recoil
		}
		outFile << std::endl;
	}
	
	outFile.close();
	if(argc >= 2){ std::cout << " Wrote output file " << argv[1] << std::endl; }
	else{ std::cout << " Wrote output file Kindist.out\n\n"; }
	
	delete[] RecoilEx;
	
	return 0;
}
