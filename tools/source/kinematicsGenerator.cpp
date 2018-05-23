
#include <iostream>

#include "vandmc_core.hpp"
#include "comConverter.hpp"

#include "TFile.h"
#include "TTree.h"

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates one output root file named 'mcarlo.root'
// fwhm_ (m) allows the use of a gaussian particle "source". If fwhm_ == 0.0, a point source is used
// angle_ (rad) allows the rotation of the particle source about the y-axis
bool GenerateKinematicsFile(const char *fname, unsigned int num_trials, comConverter *conv, const std::string &title="Kinematics File"){
	if(!conv){ return false; }
	double labAngle;
	double phiAngle;
	double comAngle;
	
	unsigned int num_trials_chunk = num_trials/10;
	unsigned int chunk_num = 1;

	TFile *file = new TFile(fname, "RECREATE");
	if(!file->IsOpen()){
		std::cout << " Error! Failed to load input file \"" << fname << "\".\n";
		return false;
	}
	
	TTree *tree = new TTree("data", title.c_str());
	tree->Branch("com", &comAngle);
	tree->Branch("lab", &labAngle);
	tree->Branch("phi", &phiAngle);

	for(unsigned int i = 0; i < num_trials; i++){
		if(i != 0 && i == num_trials_chunk*chunk_num){ // Print a status update.
			std::cout << "  " << (chunk_num++)*10 << "% - " << "Generated " << i << " of " << num_trials << " total events (" << i*100.0/num_trials << "%)\n";
		}
	
		// Generate a uniformly distributed random point on the unit sphere in the center-of-mass frame.
		UnitSphereRandom(comAngle, phiAngle); // In the CM frame.
		labAngle = conv->convertEject2lab(comAngle);
		
		tree->Fill();
	}

	file->cd();
	tree->Write();
	file->Close();
	delete file;
	
	return true;
}

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " <relfile> <ofname> [numEvents=1E6] [title=\"Kinematics File\"]\n";
}

int main(int argc, char *argv[]){
	if(argc < 3){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 2, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	comConverter *conv = new comConverter(argv[1]);
	
	unsigned int Nwanted = 1E6;
	if(argc >= 4){
		Nwanted = strtoul(argv[3], NULL, 0);
	}
	
	std::string title = "Kinematics File";
	if(argc >= 5){
		title = std::string(argv[4]);
	}	

	std::cout << " Generating " << Nwanted << " Monte Carlo events...\n";
	GenerateKinematicsFile(argv[2], Nwanted, conv, title);

	std::cout << " Finished generating kinematics Monte Carlo file...\n";
	
	return 0;
}
