
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

const double Rmin = 0.011; // m
const double Rmax = 0.035; // m
const double pi = 3.1415926536; // rad

const int numStrips = 24;
const int numSectors = 8;
const double stripWidth = (Rmax-Rmin)/numStrips; // m
const double sectorAngle = 2*pi/numSectors; // rad

/// Return the local polar coordinates of a hit on the surface of the detector.
bool getLocalPosition(const double &x_, const double &y_, double &r, double &theta){

	// Compute the radius on the local detector face.
	r = std::sqrt(x_*x_ + y_*y_);

	// Check that the radius of the event is within the active detector area.
	if(r < Rmin || r > Rmax) return false;

	// Compute the polar angle on the local detector face.
	theta = pi + std::atan2(y_, x_);

	return true;
}

/// Indexed from zero.
unsigned short getStrip(const double &r_){
	return (unsigned short)((r_-Rmin)/stripWidth);
}

/// Indexed from zero.
unsigned short getSector(const double &theta_){
	return (unsigned short)(theta_/sectorAngle);
}

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [mcarloFile]\n";
}

int main(int argc, char *argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	TFile *infile = new TFile(argv[1], "READ");

	if(!infile || !infile->IsOpen()){
		std::cout << " Error! Failed to open input file \"" << argv[1] << "\"!\n";
		return 1;
	}

	TTree *intree = (TTree*)infile->Get("data");

	if(!intree){
		std::cout << " Error! Failed to load input tree from file!\n";
		infile->Close();
		delete infile;
		return 2;
	}

	std::vector<double> hitX, hitY;
	TBranch *hitX_b = NULL;
	TBranch *hitY_b = NULL;

	intree->SetMakeClass(1);
	intree->SetBranchAddress("face1_hitX", &hitX, &hitX_b);
	intree->SetBranchAddress("face1_hitY", &hitY, &hitY_b);

	if(!hitX_b || !hitY_b){
		std::cout << " Error! Failed to load branches from input tree!\n";
		infile->Close();
		delete infile;
		return 3;
	}

	TFile *outfile = new TFile("strips.root", "RECREATE");
	TTree *outtree = new TTree("data", "Strips tree");

	double x, y;
	double r, theta;
	unsigned short sector, strip;

	outtree->Branch("x", &x);
	outtree->Branch("y", &y);
	outtree->Branch("r", &r);
	outtree->Branch("theta", &theta);
	outtree->Branch("sector", &sector);
	outtree->Branch("strip", &strip);

	for(int entry = 0; entry < intree->GetEntries(); entry++){
		intree->GetEntry(entry);
		for(size_t i = 0; i < hitX.size(); i++){
			x = hitX.at(i);
			y = hitY.at(i);

			// Convert to coordinates in the local frame of the face of the detector.
			if(!getLocalPosition(x, y, r, theta)) continue;

			sector = getSector(theta);
			strip = getStrip(r);

			outtree->Fill();
		}
	}

	outfile->cd();
	outtree->Write();
	outfile->Close();

	infile->Close();

	delete outfile;
	delete infile;

	return 0;
}
