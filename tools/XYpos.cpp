// XYpos.cpp
// Cory Thornsberry
// Oct. 27th, 2014
// Calculate the XY position spectrum for a VIKAR data set

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cmath>

#include "vikar_core.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"

#define UNCX 0.001 // 1 mm resolution on the x-grid
#define UNCY 0.001 // 1 mm resolution on the y-grid

#define VANDLE_UNCT 0.03092 // 0.03092 rad angular resolution for a VANDLE bar at 0.5 m from the target
#define VANDLE_UNCY 0.1 // 10 cm resolution for vertical position inside a VANDLE bar

struct VBar{
	unsigned int ID;
	double R, theta, phi;
	
	VBar(unsigned int ID_, double R_, double theta_, double phi_){
		ID = ID_; R = R_; theta = theta_; phi = phi_;
	}
};

bool GetTheta(const std::vector<VBar> &bars_, unsigned int ID_, double &theta){
	for(std::vector<VBar>::const_iterator iter = bars_.begin(); iter != bars_.end(); iter++){
		if(ID_ == iter->ID){ 
			theta = iter->theta;
			return true;
		}
	}
	return false;
}

int main(int argc, char *argv[]){
	gSystem->Load("libTree");
	if(argc < 6){ 
		std::cout << " Error! Invalid number of arguments\n";
		std::cout << " Syntax: ./XYpos filename treename Xdet Ydet bar_data {debug}\n";
		return 1; 
	}
	
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

	bool debug = false;
	if(strcmp(argv[6], "debug") == 0){
		debug = true;
		std::cout << " Debugging...\n";
	}
	
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
	
	std::vector<unsigned int> eloc, rloc;
	std::vector<unsigned int>::iterator iterL;
	
	std::vector<double> faceX, faceY, hitY;
	std::vector<double>::iterator iterX, iterY;

	// Ejectile Branches
	TBranch *hitY_b, *eloc_b;	
	tree->SetBranchAddress("eject_faceY", &hitY, &hitY_b);
	if(!hitY_b){ 
		std::cout << " Error! Failed to load branch 'eject_hitY'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress("eject_loc", &eloc, &eloc_b);
	if(!eloc_b){ 
		std::cout << " Error! Failed to load branch 'eject_loc'\n"; 
		file->Close();
		return 1;
	}

	// Recoil Branches
	TBranch *faceX_b, *faceY_b, *rloc_b;	
	tree->SetBranchAddress("recoil_faceX", &faceX, &faceX_b);
	if(!faceX_b){ 
		std::cout << " Error! Failed to load branch 'recoil_faceX'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress("recoil_faceY", &faceY, &faceY_b);
	if(!faceY_b){ 
		std::cout << " Error! Failed to load branch 'recoil_faceY'\n"; 
		file->Close();
		return 1;
	}
	tree->SetBranchAddress("recoil_loc", &rloc, &rloc_b);
	if(!rloc_b){ 
		std::cout << " Error! Failed to load branch 'recoil_loc'\n"; 
		file->Close();
		return 1;
	}

	unsigned int Xdet = (unsigned int)atoi(argv[3]);
	unsigned int Ydet = (unsigned int)atoi(argv[4]);
	std::cout << " X Grid Location: " << Xdet << std::endl;
	std::cout << " Y Grid Location: " << Ydet << std::endl;
	
	std::ofstream output;
	if(debug){ output.open("XYpos.out"); }
	
	std::ifstream bar_file(argv[5]);
	if(!bar_file.good()){
		std::cout << " Error! Failed to load detector input file '" << argv[6] << "'\n";
		file->Close();
		bar_file.close();
		return 1;
	}
		
	std::vector<VBar> detectors;
	unsigned int bar_id;
	double bar_R, bar_theta, bar_phi;
	while(true){
		bar_file >> bar_id >> bar_R >> bar_theta >> bar_phi;
		if(bar_file.eof()){ break; }
		detectors.push_back(VBar(bar_id, bar_R, bar_theta, bar_phi));
	}
	
	TH2F *hist = new TH2F("hist", "Ion Chamber X vs. Y", 254, -0.0381, 0.0381, 254, -0.0381, 0.0381);
	hist->GetXaxis()->SetTitle("X (0.3 mm/bin)");
	hist->GetYaxis()->SetTitle("Y (0.3 mm/bin)");
	hist->SetStats(false);
	
	Vector3 temp_vectors[4];
	Vector3 ejectP[4];
	Vector3 recoilP[4];
	double X, Y;
	
	double eject_Y;
	double temp_theta;
	unsigned int eject_loc;
	
	Matrix3 Vandle_Left = Matrix3(-VANDLE_UNCT, 0.0);
	Matrix3 Vandle_Right = Matrix3(VANDLE_UNCT, 0.0);
	
	unsigned int event_counts = 0;
	unsigned int num_entries = tree->GetEntries();
	if(debug && 10 < tree->GetEntries()){
		num_entries = 10;
	}
	
	std::cout << " Processing " << num_entries << " tree entries\n";
	for(unsigned int i = 0; i < num_entries; i++){
		tree->GetEntry(i);
		X = -1.0;
		Y = -1.0;
		
		// Ejectile information (1 per event)
		eject_Y = (*hitY.begin());
		eject_loc = (*eloc.begin());
		
		if(!GetTheta(detectors, eject_loc, temp_theta)){
			std::cout << " Problem finding ejectile detector #" << eject_loc << std::endl;
			continue;
		}

		if(debug){
			double temp_phi = std::atan(eject_Y/0.5);
			std::cout << i << ": eject_Y = " << eject_Y << ",  eject_loc = " << eject_loc;
			std::cout << ",  eject_theta = " << temp_theta*rad2deg << ",  eject_phi = " << temp_phi*rad2deg << std::endl;
		}
				
		// Due to the VANDLE resolution, the intersect can occur anywhere inside a square with vertices at...
		ejectP[0] = Vector3(0.5, temp_theta, 0.0);
		ejectP[1] = Vector3(0.5, temp_theta, 0.0);
		ejectP[2] = Vector3(0.5, temp_theta, 0.0);
		ejectP[3] = Vector3(0.5, temp_theta, 0.0);
		
		for(unsigned int j = 0; j < 4; j++){ Sphere2Cart(ejectP[j], temp_vectors[j]); }		
		
		Vandle_Left.Transform(temp_vectors[0]);
		Vandle_Left.Transform(temp_vectors[2]);
		Vandle_Right.Transform(temp_vectors[1]);
		Vandle_Right.Transform(temp_vectors[3]);
		
		temp_vectors[0].axis[1] = eject_Y + VANDLE_UNCY/2.0;
		temp_vectors[1].axis[1] = eject_Y + VANDLE_UNCY/2.0;
		temp_vectors[2].axis[1] = eject_Y -VANDLE_UNCY/2.0;
		temp_vectors[3].axis[1] = eject_Y -VANDLE_UNCY/2.0;

		// Radii are approximate, we shouldn't need them anyway
		for(unsigned int j = 0; j < 4; j++){ 
			Cart2Sphere(temp_vectors[j], ejectP[j]); 
			if(debug){ output << temp_vectors[j].axis[0] << "\t" << temp_vectors[j].axis[1] << "\t" << temp_vectors[j].axis[2] << "\n"; }
		}	
		
		// Recoil information (multiple per event)
		for(iterX = faceX.begin(), iterY = faceY.begin(), iterL = rloc.begin(); 
		  iterX != faceX.end() && iterY != faceY.end() && iterL != rloc.end(); iterX++, iterY++, iterL++){
		  	if(*iterL == Xdet){ X = *iterX; }
		  	else if(*iterL == Ydet){ Y = *iterY; }
		}
		
		if(X != -1.0 && Y != -1.0){
			// Assume the grids are sufficiently close together such that the recoil cone from the
			// target is essentially a cylinder between the X and Y grids of the ion chamber
			
			// Due to the resolution on the IC grids, the intersect can occur anywhere inside a square with vertices at...
			temp_vectors[0] = Vector3(X-UNCX, Y+UNCY, 0.105); // Top left of the bounding box
			temp_vectors[1] = Vector3(X+UNCX, Y+UNCY, 0.105); // Top right of the bounding box
			temp_vectors[2] = Vector3(X+UNCX, Y-UNCY, 0.105); // Bottom right of the bounding box
			temp_vectors[3] = Vector3(X-UNCX, Y-UNCY, 0.105); // Bottom left of the bounding box
			
			// Convert the ion-chamber intersect to spherical to extract the reaction information
			for(unsigned int j = 0; j < 4; j++){ 
				Cart2Sphere(temp_vectors[j], recoilP[j]); 
				if(debug){ output << temp_vectors[j].axis[0] << "\t" << temp_vectors[j].axis[1] << "\t" << temp_vectors[j].axis[2] << "\n"; }
			}

			
			
			hist->Fill(X, Y);
			event_counts++;
		}
	}
	std::cout << "  Done! Found " << event_counts << " events\n";
	
	/*TCanvas *can = new TCanvas("can");
	can->cd();
	hist->Draw("COLZ");
	can->WaitPrimitive();
	can->Close();*/
	
	file->Close();
	rootapp->Delete();
	
	if(debug){ output.close(); }
	
	return 0;
}
