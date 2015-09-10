// mc_viewer.cpp
// C. Thornsberry
// Sept. 10th, 2015

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TApplication.h"

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [filename] <options>\n";
	std::cout << "   Available options:\n";
	std::cout << "    --box | Draw using the 'BOX' draw option.\n";
	std::cout << "    --iso | Draw using the 'ISO' draw option.\n";
}

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	int index = 2;
	std::string draw_opt = "";
	while(index < argc){
		if(strcmp(argv[index], "--box") == 0){ draw_opt = "BOX"; }
		else if(strcmp(argv[index], "--iso") == 0){ draw_opt = "ISO"; }
		else{ 
			std::cout << " Error! Unrecognized option '" << argv[index] << "'!\n";
			help(argv[0]);
			return 1;
		}
		index++;
	}

	// Variables for root graphics
	TApplication* rootapp = new TApplication("rootapp",0,NULL);
	gSystem->Load("libTree");
	
	TFile *file = new TFile(argv[1], "READ");
	if(!file->IsOpen()){
		std::cout << " Failed to load the input file '" << argv[1] << "'\n";
		return 1;
	}
	TTree *tree = (TTree*)file->Get("VIKAR");
	if(!tree){
		std::cout << " Failed to load the input tree 'Pixie16'\n";
		file->Close();
		return 1;
	}
	
	TCanvas *can = new TCanvas("can", "MC Viewer");
	can->cd();
	
	tree->Draw("face1_hitY:face1_hitZ:face1_hitX","",draw_opt.c_str());

	can->WaitPrimitive();
	
	can->Close();
	file->Close();
	
	rootapp->Delete();
	
	return 0;
}
