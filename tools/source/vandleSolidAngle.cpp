#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "vandmc_core.h"
#include "detectors.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TNamed.h"

// List the ratio of each bin in two 1-d histograms.
void binRatio(TH1 *h1_, TH1 *h2_, TH1 *h3_, TH1 *h4_, const std::vector<double> &comAngles_, std::ostream &out=std::cout){
	out << "\nbin\tlabCenter\tcomLow\tcomCenter\tcomErr\tN1\tN2\tN3\tN4\tbinSolidAngle\tSA0\tSA1\tSA2\tdSAtheta(%)\n";
	out.precision(6);
	double binWidth;
	double content1, content2;
	double content3, content4;
	double angleLow, angleHigh;
	double binSolidAngle;
	double solidAngle[3];
	for(int i = 1; i <= h1_->GetNbinsX(); i++){
		content1 = h1_->GetBinContent(i);
		content2 = h2_->GetBinContent(i);
		content3 = h3_->GetBinContent(i);
		content4 = h4_->GetBinContent(i);

		angleLow = comAngles_.at(i-1);
		angleHigh = comAngles_.at(i);

		if(angleLow > angleHigh){ // Inverse kinematics.
			double tempValue = angleLow;
			angleLow = angleHigh;
			angleHigh = tempValue;
		}

		binWidth = angleHigh-angleLow;
		binSolidAngle = 2*pi*(std::cos(angleLow)-std::cos(angleHigh));
		
		solidAngle[0] = binSolidAngle*(content2/content1);
		solidAngle[1] = binSolidAngle*(content3/content1);
		solidAngle[2] = binSolidAngle*(content4/content1);
		
		double maxDeviation = std::max(std::fabs(solidAngle[1]-solidAngle[0]), std::fabs(solidAngle[1]-solidAngle[2]));
		
		out << i << "\t" << h1_->GetBinCenter(i)*rad2deg << "\t" << angleLow*rad2deg << "\t";
		out << (angleLow+binWidth/2)*rad2deg << "\t" << (binWidth/2)*rad2deg << "\t";
		out << content1 << "\t" << content2 << "\t" << content3 << "\t" << content4 << "\t";
		out << std::fixed << binSolidAngle << "\t" << solidAngle[0] << "\t";
		out << solidAngle[1] << "\t" << solidAngle[2] << "\t" << 100*maxDeviation/solidAngle[1] << std::endl;
		out.unsetf(std::ios_base::floatfield);
	}
}

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates one output root file named 'mcarlo.root'
// fwhm_ (m) allows the use of a gaussian particle "source". If fwhm_ == 0.0, a point source is used
// angle_ (rad) allows the rotation of the particle source about the y-axis
unsigned int TestDetSetup(TTree *comTree, const std::vector<double> &barAngles, const std::vector<double> &angleBins, const std::vector<double> &comBins, const double &radius_){
	if(!comTree){ return 0; }
	double labAngle;
	double phiAngle;
	double comAngle;
	unsigned int totalInAngleBin=0;
	unsigned int totalDetected=0;
	Vector3 flight_path;
	Vector3 temp_ray;
	
	TBranch *lab_b=NULL, *phi_b=NULL, *com_b=NULL;
	comTree->SetBranchAddress("lab", &labAngle, &lab_b);
	comTree->SetBranchAddress("phi", &phiAngle, &phi_b);
	comTree->SetBranchAddress("com", &comAngle, &com_b);
	
	if(!lab_b || !phi_b || !com_b){
		return 0;
	}
	
	TH1F *h1 = new TH1F("h1", "Ungated", angleBins.size()-1, angleBins.data());
	TH1F *h2 = new TH1F("h2", "Detector Gated (mask 0)", angleBins.size()-1, angleBins.data());
	TH1F *h3 = new TH1F("h3", "Detector Gated (mask 1)", angleBins.size()-1, angleBins.data());
	TH1F *h4 = new TH1F("h4", "Detector Gated (mask 2)", angleBins.size()-1, angleBins.data());
	
	unsigned int Nentries = comTree->GetEntries();
	unsigned int num_trials_chunk = comTree->GetEntries()/10;
	unsigned int chunk_num = 1;
	double t, x, y;

	const double width = 0.03; // m
	const double barHalfAngle = std::asin(width/(2*radius_));

	// Compute mask angles.
	const double dTheta = 2; // deg
	const double angles[3] = {barHalfAngle-dTheta*deg2rad, barHalfAngle, barHalfAngle+dTheta*deg2rad};

	for(unsigned int i=0; i < Nentries; i++){
		if(i != 0 && i == num_trials_chunk*chunk_num){ // Print a status update.
			std::cout << "  " << (chunk_num++)*10 << "% - " << "Detected " << i << " of " << Nentries << " total events (" << i*100.0/Nentries << "%)\n";
		}
	
		// Get an event from the tree.
		comTree->GetEntry(i);
	
		// Fill the ungated histogram.
		h1->Fill(labAngle);
	
		// Check CM angle. To speed things up.
		if(comAngle >= pi/2) continue;

		// Convert from spherical to cartesian.
		Sphere2Cart(1.0, labAngle, phiAngle, temp_ray);

		totalInAngleBin++;

		t = radius_/std::sqrt(temp_ray.axis[0]*temp_ray.axis[0]+temp_ray.axis[2]*temp_ray.axis[2]);
		x = t*temp_ray.axis[0];
		y = t*temp_ray.axis[1];
		//z = t*temp_ray.axis[2];
		
		// Check for events which do not intersect VANDLE.
		if(y >= -0.3 && y <= 0.3 && x >= 0){
			// Check for detector hit.
			int loc = -1;
			int mask = -1;
			for(size_t i = 0; i < barAngles.size(); i++){
				if(labAngle >= barAngles[i]-angles[2] && labAngle <= barAngles[i]+angles[2]){
					loc = i;
					if(labAngle < barAngles[i]){
						if(labAngle < barAngles[i]-angles[1]) mask = 0;
						else if(labAngle < barAngles[i]-angles[0]) mask = 1;
						else mask = 2;
					}
					else{
						if(labAngle >= barAngles[i]+angles[1]) mask = 4;
						else if(labAngle >= barAngles[i]+angles[0]) mask = 3;
						else mask = 2;
					}
					break;
				}
			}

			if(loc >= 0 && mask >= 0){
				if(mask <= 2){ // Bar mask -dTheta
					h2->Fill(labAngle);
				}
				if(mask >=1 && mask <= 3){ // Bar mask
					h3->Fill(labAngle);
					totalDetected++;
				}
				if(mask >=2 && mask <= 4){ // Bar mask +dTheta
					h4->Fill(labAngle);
				}
			}
		}
	}

	// Compute the solid angle per bin.
	binRatio(h1, h2, h3, h4, comBins); 

	delete h1;
	delete h2;
	delete h3;
	delete h4;

	return totalDetected;
}

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " <detfile> <binfile> [MCfile]\n";
}

int main(int argc, char *argv[]){
	if(argc < 3){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 2, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	std::vector<Primitive*> detectors;

	std::cout << " Reading in NewVIKAR detector setup file...\n";
	int Ndet = ReadDetFile(argv[1], detectors);
	if(Ndet < 0){
		std::cout << " Error: failed to load detector setup file!\n";
		return 1; 
	}
	else if(Ndet == 0){ std::cout << " Error: Found no detectors in the detector setup file!\n"; }

	std::cout << "  Loaded " << Ndet << " detectors from file.\n";

	std::cout << " Reading angular bin file...\n";
	std::ifstream binFile(argv[2]);
	if(!binFile.good()){
		std::cout << " Error: Failed to load angular bin file!\n";
		return 1;
	}
	
	double labTheta, comTheta;
	std::vector<double> labBins;
	std::vector<double> comBins;
	while(true){
		binFile >> labTheta >> comTheta;
		if(binFile.eof()) break;
		labBins.push_back(labTheta*deg2rad);
		comBins.push_back(comTheta*deg2rad);
	}
	binFile.close();

	std::cout << "  Loaded " << labBins.size()-1 << " angular bins.\n"; 

	std::vector<double> detectorAngles;
	for(std::vector<Primitive*>::iterator iter = detectors.begin(); iter != detectors.end(); iter++){
		detectorAngles.push_back((*iter)->GetTheta());
	}

	double radius;
	std::cout << " Enter detector radius (m): "; std::cin >> radius;
	
	TFile *file = NULL;
	if(argc >= 4){
		std::string rxnFilename(argv[3]);
		file = new TFile(argv[3], "READ");
		if(!file->IsOpen()){
			std::cout << " Error: failed to load input kinematics file \"" << argv[3] << "\"!\n";
			return 1;
		}
	}
	
	TTree *tree = NULL;
	tree = (TTree*)file->Get("data");
	
	if(!file->IsOpen()){
		std::cout << " Error: failed to load input kinematics file \"" << argv[3] << "\"!\n";
		return 1;
	}	

	std::cout << std::endl;

	std::cout << "  Performing Monte Carlo test on VANDLE detectors...\n";
	unsigned int totalDetected = TestDetSetup(tree, detectorAngles, labBins, comBins, radius);

	std::cout << "  Detected " << totalDetected << " ejectile events in " << tree->GetEntries() << " trials (" << 100.0*totalDetected/tree->GetEntries() << "%)\n\n";

	std::cout << " Finished geometric efficiency test on VANDLE setup...\n";

	return 0;
}
