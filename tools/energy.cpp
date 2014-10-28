// energy.cpp
// Cory Thornsberry
// Oct. 16th, 2014
// Calculate the energy spectrum for a VIKAR data set

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cmath>

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"

#define MASS 939.565378 // MeV/c^2
#define PI 3.1415926540

void Cart2Sphere(double x, double y, double z, double &r, double &theta, double &phi){ 
	if(x == 0.0 && y == 0.0 && z == 0.0){
		r = 0.0; theta = 0.0; phi = 0.0;
	}
	else{
		r = std::sqrt(x*x + y*y + z*z);
		theta = std::acos(z/r);
		
		if(x == 0.0 && y == 0.0){ phi = 0.0; }
		else{ 
			double temp = std::sqrt(x*x + y*y); 
			if(x >= 0.0){ phi = std::acos(y/temp); }
			else{ phi = 2.0*PI - std::acos(y/temp); }
		}
	}
}

int main(int argc, char *argv[]){
	gSystem->Load("libTree");
	if(argc < 4){ 
		std::cout << " Error! Invalid number of arguments\n";
		std::cout << " Syntax: ./energy filename treename branch_prefix\n";
		return 1; 
	}
	
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	
	TFile *file = new TFile(argv[1], "READ");
	if(!file->IsOpen()){
		std::cout << " Error! Failed to load the input file '" << argv[1] << "'\n";
		return 1;
	}
	
	TTree *tree = (TTree*)file->Get(argv[2]);
	if(!tree){
		std::cout << " Error! Failed to load the input tree '" << argv[2] << "'\n";
		return 1;
	}
	tree->SetMakeClass(1);
			
	std::vector<double> tof, hitX, hitY, hitZ;
	std::vector<double>::iterator iterT, iterX, iterY, iterZ;
	TBranch *tof_b, *hitX_b, *hitY_b, *hitZ_b;
	
	std::stringstream stream; stream << argv[3];
	std::string tof_branchname = stream.str()+"_tof";
	std::string hitX_branchname = stream.str()+"_hitX";
	std::string hitY_branchname = stream.str()+"_hitY";
	std::string hitZ_branchname = stream.str()+"_hitZ";
	
	tree->SetBranchAddress(tof_branchname.c_str(), &tof, &tof_b);
	if(!tof_b){ 
		std::cout << " Warning! Failed to load branch '" << tof_branchname << "'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress(hitX_branchname.c_str(), &hitX, &hitX_b);
	if(!hitX_b){ 
		std::cout << " Warning! Failed to load branch '" << hitX_branchname << "'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress(hitY_branchname.c_str(), &hitY, &hitY_b);
	if(!hitY_b){ 
		std::cout << " Warning! Failed to load branch '" << hitY_branchname << "'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress(hitZ_branchname.c_str(), &hitZ, &hitZ_b);
	if(!hitZ_b){ 
		std::cout << " Warning! Failed to load branch '" << hitZ_branchname << "'\n"; 
		file->Close();
		return 1;
	}
	
	TH2F *hist = new TH2F("hist", "Neutron E vs. Lab Angle", 180, 0, 180, 100, 0, 12);
	hist->GetXaxis()->SetTitle("Lab Angle (1 deg/bin)");
	hist->GetYaxis()->SetTitle("Neutron Energy (0.12 MeV/bin)");
	hist->SetStats(false);
	
	double r, theta, phi, energy, time;
	unsigned int event_counts = 0;
	std::cout << " Processing " << tree->GetEntries() << " tree entries\n";
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		for(iterT = tof.begin(), iterX = hitX.begin(), iterY = hitY.begin(), iterZ = hitZ.begin();
		  iterT != tof.end() && iterX != hitX.end() && iterY != hitY.end() && iterZ != hitZ.end();
		  iterT++, iterX++, iterY++, iterZ++){
			Cart2Sphere(*iterX, *iterY, *iterZ, r, theta, phi);
			time = *iterT;
			energy = 0.5*(MASS/(9E-2))*r*r/(time*time); // Particle energy in MeV
			hist->Fill(theta*180.0/PI, energy);
			event_counts++;
		}
	}
	std::cout << "  Done! Found " << event_counts << " events\n";
	
	TCanvas *can = new TCanvas("can");
	can->cd();
	hist->Draw("COLZ");
	can->WaitPrimitive();
	can->Close();
	
	file->Close();
	rootapp->Delete();
	
	return 0;
}
