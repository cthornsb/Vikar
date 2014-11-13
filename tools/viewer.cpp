// PulseViewer.cpp
// C. Thornsberry
// Aug. 25th, 2014
// Load detector pulses and view them in sequence
// SYNTAX: Viewer {root_filename} {root_treename} {root_branchname} {skip#}

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TApplication.h"

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>

double PIXIE_TIME_RES = 4.0; // In ns

// For compilation
int main(int argc, char* argv[]){
	if(argc < 5){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: Viewer {filename} {treename} {branch} {skip#}\n";
		return 1;
	}
	
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	unsigned int skip = atol(argv[4]);
	std::cout << " Showing every " << skip << " pulses\n";

	// Branch variables
	std::vector<int> wave;
	std::vector<double> energy;
	unsigned int wave_size = 0;
		
	TFile *file = new TFile(argv[1], "READ");
	if(!file->IsOpen()){
		std::cout << " Failed to load the input file '" << argv[1] << "'\n";
		return 1;
	}
	TTree *tree = (TTree*)file->Get(argv[2]);
	if(!tree){
		std::cout << " Failed to load the input tree '" << argv[2] << "'\n";
		file->Close();
		return 1;
	}
	tree->SetMakeClass(1);

	TBranch *b_wave;
	tree->SetBranchAddress(argv[3], &wave, &b_wave);
	
	if(!b_wave){
		std::cout << " Failed to load the input branch '" << argv[3] << "'\n";
		file->Close();
		return 1;
	}

	// Get the pulse size
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if(wave.size() == 0){ continue; }
		else{ 
			wave_size = wave.size();
			break; 
		}
	}
	std::cout << " Using wave size " << wave_size << std::endl;

	TCanvas *can = new TCanvas("can", "canvas");
	can->cd();

	// Variables for TGraph
	double *x = new double[wave_size];
	const double *x_val = new double[wave_size];
	const double *y_val = new double[wave_size];
	
	TGraph *graph = new TGraph(wave_size, x_val, y_val);
	for(unsigned int i = 0; i < wave_size; i++){
		x[i] = i*PIXIE_TIME_RES;
	}

	unsigned int index = 0;
	std::cout << " Processing " << tree->GetEntries() << " entries\n";
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if(i % (skip*10) == 0){ std::cout << " Entry no. " << i << std::endl; }
		if(i % skip == 0){ 
			index = 0;
			for(std::vector<int>::iterator iter = wave.begin(); iter != wave.end(); iter++){
				if(index >= wave_size){ break; } // Reached maximum graph container size
				graph->SetPoint(index, x[index], (*iter));
				index++;
			}
			graph->Draw();
			can->Update();
			sleep(1);
		}
	}

	can->Close();
	file->Close();
	graph->Delete();
	
	rootapp->Delete();
	
	return 0;
}

// For CINT
int Viewer(int argc, char* argv[]){ return main(argc, argv); }
