// dump.cpp
// Cory Thornsberry
// Oct. 16th, 2014
// Calculate the energy spectrum for a VIKAR data set

#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

int main(int argc, char *argv[]){
	gSystem->Load("libTree");
	if(argc < 4){ 
		std::cout << " Error! Invalid number of arguments\n";
		std::cout << " Syntax: ./dump filename treename N branch1 branch2... branchN\n";
		return 1; 
	}

	unsigned int num_vars = (unsigned int)atoi(argv[3]);
	if(num_vars == 0){
		std::cout << " Error! Attempt to load zero branches from input file\n";
		return 1;
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
		
	if((unsigned int)argc < (4+num_vars)){
		std::cout << " Error! Received " << argc-4 << " branch names, expected " << num_vars << "\n";
		file->Close();
		return 1;
	}
	
	std::vector<double> vars[num_vars];
	std::vector<double>::iterator iters[num_vars];
	TBranch *branches[num_vars];
	unsigned int good_branch = 0;
	bool switches[num_vars];
	bool temp = false;
	
	std::cout << " Loading " << num_vars << " root branches\n";
	for(unsigned int i = 0; i < num_vars; i++){
		tree->SetBranchAddress(argv[4+i], &vars[i], &branches[i]);
		if(!branches[i]){ 
			std::cout << " Warning! Failed to load branch '" << argv[4+i] << "'\n"; 
			switches[i] = false;
		}
		else{ 
			good_branch++;
			switches[i] = true; 
			if(!temp){ temp = true; }
		}
	}
	
	if(!temp){
		std::cout << " Error! Failed to load any input branches\n";
		file->Close();
		return 1;
	}
	
	std::cout << " Successfully loaded " << good_branch << " branches\n";
	
	std::ofstream output_file("dump.out");
	if(!output_file.good()){
		std::cout << " Error! Failed to open output file\n";
		file->Close();
		return 1;
	}
	
	bool good_entry;
	unsigned int size, count;
	unsigned int event_counts = 0;
	std::cout << " Processing " << tree->GetEntries() << " tree entries\n";
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		size = vars[0].size();
		
		good_entry = true;
		for(unsigned int j = 0; j < num_vars; j++){
			if(!switches[j]){ continue; }
			if(vars[j].size() != size){ 
				std::cout << "  vars[" << j << "].size() != " << size << std::endl;
				good_entry = false; 
				break;
			}
			iters[j] = vars[j].begin();
		}
		
		if(good_entry){
			count = 0;
			while(count < size){
				for(unsigned int j = 0; j < num_vars; j++){
					if(switches[j]){ output_file << *iters[j]; }
					else{ output_file << "-1"; }
					if(j != num_vars-1){ output_file << "\t"; }
					iters[j]++;
				}
				output_file << "\n";
				event_counts++;
				count++;
			}
		}
	}
	std::cout << "  Done! Wrote " << event_counts << " events to file\n";
	
	output_file.close();
	file->Close();
	
	return 0;
}
