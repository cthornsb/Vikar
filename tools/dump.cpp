// dump.cpp
// Cory Thornsberry
// Oct. 16th, 2014
// Dump entries from a root file with branches of vectors to a text file

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
	
	unsigned int num_vars = (unsigned int)atoi(argv[3]);
	std::cout << " Loading " << num_vars << " root branches\n";
	
	if((unsigned int)argc < (4+num_vars)){
		std::cout << " Error! Received " << argc-4 << " branch names, expected " << 4+num_vars << "\n";
		file->Close();
		return 1;
	}
	
	std::vector<double> vars[num_vars];
	std::vector<double>::iterator iters[num_vars];
	TBranch *branches[num_vars];
	
	for(unsigned int i = 0; i < num_vars; i++){
		tree->SetBranchAddress(argv[4+i], &vars[i], &branches[i]);
		if(!branches[i]){ std::cout << " Warning! Failed to load branch '" << argv[4+i] << "'\n"; }
	}
	
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
			if(vars[j].size() != size){ 
				std::cout << "  Problem with entry no. " << i << std::endl;
				good_entry = false; 
				break;
			}
			iters[j] = vars[j].begin();
		}
		
		if(good_entry){
			count = 0;
			while(count < size){
				output_file << *iters[0];
				for(unsigned int j = 1; j < num_vars; j++){
					output_file << "\t" << *iters[j];
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
