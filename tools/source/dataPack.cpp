#include "dataPack.hpp"

#include "TFile.h"
#include "TTree.h"

dataPack::dataPack(){
	file = NULL;
	tree = NULL;
	init = false;
}

dataPack::dataPack(std::string fname_, unsigned int Nbins_, bool write_rxn_/*=false*/){
	Open(fname_, write_rxn_);
}

dataPack::~dataPack(){
	if(init){ Close(); }
}

bool dataPack::IsInit(){ return init; }

bool dataPack::Open(std::string fname_, bool write_rxn_){
	if(init){ return false; }

	file = new TFile(fname_.c_str(), "RECREATE");
	
	if(!file->IsOpen()){
		init = false;
		delete file;
	}
	
	tree = new TTree("data", "Monte carlo detector efficiency tree");
	
	tree->Branch("mcarlo", &MCARLOdata);
	if(write_rxn_){ 
		tree->Branch("labTheta", &labTheta);
		tree->Branch("labPhi", &labPhi);
		tree->Branch("comAngle", &comAngle);
		tree->Branch("labBin", &labBin); 
	}
	
	return (init = true);
}

bool dataPack::Close(){
	if(!init){ return false; }
	
	file->cd();
	tree->Write();
	file->Close();
	
	init = false;
	
	return true;
}

void dataPack::Zero(){
	MCARLOdata.Zero();
}
