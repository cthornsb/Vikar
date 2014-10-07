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
	double *RecoilEx, *EjectTheta;
	double *RecoilTheta, *Eeject, *Erecoil;
	double *RecoilTheta2, *Eeject2, *Erecoil2;
	/*std::cout << "  Enter beam Energy (MeV): "; std::cin >> Ebeam;
	std::cout << "  Enter beam Mass (A): "; std::cin >> Mbeam;
	std::cout << "  Enter target Mass (A): "; std::cin >> Mtarg;
	std::cout << "  Enter g.s. Q-value (MeV): "; std::cin >> Q;
	std::cout << "  Enter ejectile Mass (A): "; std::cin >> Meject;
	std::cout << "  Enter number of states: "; std::cin >> num_states;
	if(num_states > 0){
		RecoilEx = new double[num_states];
		EjectTheta = new double[num_states];
		RecoilTheta = new double[num_states];
		for(unsigned int i = 0; i < num_states; i++){
			std::cout << "   Enter recoil excitation for state " << i+1 << " (MeV): "; 
			std::cin >> RecoilEx[i];
			EjectTheta[i] = 0.0;
			RecoilTheta[i] = 0.0;
		}
	}
	else{
		RecoilEx = new double[1]; RecoilEx[0] = 0.0;
		EjectTheta = new double[1]; EjectTheta[0] = 0.0;
		RecoilTheta = new double[1]; RecoilTheta[0] = 0.0;
	}

	std::cout << "  Enter start Lab angle (deg): "; std::cin >> start_angle; start_angle *= deg2rad;
	std::cout << "  Enter stop Lab angle (deg): "; std::cin >> stop_angle; stop_angle *= deg2rad;
	std::cout << "  Enter number of steps: "; std::cin >> num_steps;*/

	Ebeam = 30.3; Mbeam = 7; Mtarg = 2; Q = -2.085;
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
	start_angle = 0.0; stop_angle = pi; num_steps = 180;

	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q, num_states, RecoilEx, 0.0);
	
	std::ofstream outFile;
	if(argc >= 2){ 
		outFile.open(argv[1]); 
		if(!outFile.good()){
			std::cout << "  Error: Failed to open the output file, check that the path is correct\n";
			delete[] RecoilEx;
			return 1;
		}
	}
	else{ outFile.open("Kindist.out"); }
	outFile << "# Kindist kinematics output file v. 1.0\n";
	outFile << "# Ebeam = " << Ebeam << " MeV\n";
	outFile << "# Mbeam = " << Mbeam << " AMU\n";
	outFile << "# Mtarg = " << Mtarg << " AMU\n";
	outFile << "# Meject = " << Meject << " AMU\n";
	outFile << "# num_states = " << num_states << "\n";
	for(unsigned int i = 0; i < num_states; i++){
		outFile << "#  RecoilEx[" << i << "] = " << RecoilEx[i] << " MeV\n";
	}
	outFile << "# start_angle = " << start_angle*rad2deg << " degrees (Lab frame)\n";
	outFile << "# stop_angle = " << stop_angle*rad2deg << " degrees (Lab frame)\n";
	outFile << "# num_steps = " << num_steps << "\n";
	outFile << "########################################################################################################################\n";
	outFile << "# CoMAngle		EjectTheta		EjectE1		EjectE2		RecoilTheta1		RecoilE1		RecoilTheta2		RecoilE2\n";
	
	double step_size = (stop_angle - start_angle)/num_steps;
	for(unsigned int i = 0; i < num_steps; i++){
		kind.FillVars(Ebeam, (start_angle+i*step_size), EjectTheta, RecoilTheta, RecoilTheta2, Eeject, Eeject2, Erecoil, Erecoil2);
		outFile << (start_angle+i*step_size);
		for(unsigned int j = 0; j < num_states; j++){
			outFile << "\t" << EjectTheta[j] << "\t" << Eeject[j] << "\t" << Eeject2[j] << "\t"; // Ejectile
			outFile << RecoilTheta[j] << "\t" << Erecoil[j] << "\t" << RecoilTheta2[j] << "\t" << Erecoil2[j]; // Recoil
		}
		outFile << std::endl;
	}
	
	outFile.close();
	if(argc >= 2){ std::cout << " Wrote output file " << argv[1] << std::endl; }
	else{ std::cout << " Wrote output file Kindist.out\n\n"; }
	
	delete[] RecoilEx;
	delete[] EjectTheta;
	delete[] RecoilTheta;
	delete[] Eeject;
	delete[] Erecoil;
	
	return 0;
}
