// Cory Thornsberry
// 08/19/2013
// Plotter

#define PLOTTER_VERSION "04.09.2014A"

#include "TCutG.h"

#include "/home/cory/Research/ROOTlib/FSExecutable.hpp"
#include "/home/cory/Research/ROOTlib/FSPlotter.hpp"

using namespace FightingSquirrel;

int main(int argc, char* argv[]){	
	std::cout << "------------------------------------------------------------\n";
	std::cout << "-                      SpecialCutter                       -\n";
	std::cout << "------------------------------------------------------------\n";

	if(argc >= 2 && strcmp(argv[1],"version") == 0){ // Display version information
		std::cout << " SpecialPlotter v" << PLOTTER_VERSION << " by Cory Thornsberry\n";
		std::cout << "  FightingSquirrel Library Dependencies\n";
		std::cout << "   Framework\tv" << FSVersions::Framework() << std::endl;
		std::cout << "   Plotter2D\tv" << FSVersions::Plotter2D() << std::endl;
		return 0;
	}

	char* chicken[0]; TApplication* rootapp = new TApplication("rootapp",0,chicken);
	std::string files[10] = {"Thresh000_main.root", "Thresh005_main.root", "Thresh010_main.root",
				 "Thresh050_main.root", "Thresh100_main.root", "Thresh200_main.root",
				 "Thresh300_main.root", "Thresh400_main.root", "Thresh500_main.root",
				 "Thresh1000_main.root"};

	// Main
	Plotter2D Plot;

	float LabTheta;
	float Nenergy;
	unsigned int counts;	

	TTree *plot_tree;
	TFile *root_files[10];
	TTree *root_trees[10];
	for(unsigned short i = 0; i < 10; i++){
		root_files[i] = new TFile(files[i].c_str(),"READ");
		if(root_files[i]->IsZombie()){
			std::cout << " Bad file: " << files[i] << std::endl;
			return 1;
		}
		
		root_trees[i] = (TTree*)root_files[i]->Get("Tree");
		if(!root_trees[i]){
			std::cout << " Bad tree in file " << files[i] << std::endl;
			return 1;
		}
		
		if(i == 0){ plot_tree = (TTree*)root_trees[0]->Clone("plot_tree"); }
		root_trees[i]->SetBranchAddress("Col03",&LabTheta);
		root_trees[i]->SetBranchAddress("Col06",&Nenergy);
		counts = 0;
	}

	Plot.Initialize(plot_tree->GetName(),plot_tree,400,twoTuple<double>(0.0,3.14159),400,twoTuple<double>(0.0,7E-7));
	std::cout << " Found " << Plot.FillHist("Col03","Col06") << " entries\n";
	Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (rad)");
	Plot.GetHist()->GetYaxis()->SetTitle("Neutron ToF (ns)");

	Plot.GetHist()->GetYaxis()->SetTitleOffset(1.2);
	Plot.MakeDrawable();
	Plot.Draw("COLZ",false);
	
	std::cout << " Make cut now...\n";
	TCutG *cut = (TCutG*)Plot.GetCanvas()->WaitPrimitive("CUTG");
	
	for(unsigned short i = 0; i < 10; i++){
		counts = 0;
		for(unsigned int j = 0; j < root_trees[i]->GetEntries(); j++){
			root_trees[i]->GetEntry(j);
			if(cut->IsInside(LabTheta,Nenergy)){ counts++; }
		}
		
		if(counts > 0){ std::cout << "  " << files[i] << "\t" << root_trees[i]->GetEntries() << "\t" << counts << "\t" << 100.0*counts/root_trees[i]->GetEntries() << "%\n"; }
		else{ std::cout << "  " << files[i] << "\t" << root_trees[i]->GetEntries() << "\t0\tN/A\n"; }
		root_files[i]->Close();
	}

	cut->Delete();

	return 0;
}
