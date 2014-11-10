#include <iostream>
#include <fstream>

#include "vikar_core.h"
#include "materials.h"
#include "detectors.h"
#include "structures.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

int main(){
	char* dummy[0]; 
	gSystem->Load("libTree");
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	
	double Ebeam0 = 30.3; // MeV
	double EbeamSpread = 4.0; // MeV
	double thetaMax = 13.3241*deg2rad; // radians
	double beamspot = 0.02; // m
	double timeRes = 3E-9; // s
	double detectorspot = beamspot+0.3*std::tan(thetaMax); // m
	double focalpoint = -beamspot/(2.0*std::tan(thetaMax)); // m
	
	std::vector<Planar*> detectors;
	unsigned int num_dets = ReadDetFile("./detectors/Phoswich.det", detectors);
	if(num_dets == 0){
		std::cout << " Error: Failed to load detector file\n";
		return 1;
	}
	
	std::cout << " Successfully loaded " << num_dets << " detectors from file\n";

	// All of the possible beam particles
	Particle particles[6];
	particles[0].SetParticle("Li6", 3, 6); // Lithium-6
	particles[1].SetParticle("Li7", 3, 7); // Lithium-7
	particles[2].SetParticle("Be7", 4, 7); // Beryillium-7
	particles[3].SetParticle("Be8", 4, 8); // Beryillium-8 (short lived)
	particles[4].SetParticle("B8", 5, 8); // Boron-8
	particles[5].SetParticle("B9", 5, 9); // Boron-9 (short lived)
	
	Material BC408;
	BC408.ReadMatFile("./materials/BC408.mat");
	std::cout << "  Beam energy: " << Ebeam0 << " MeV\n";
	std::cout << "  Energy spread: " << EbeamSpread << " MeV\n";
	std::cout << "  Theta max: " << thetaMax << " rad\n";
	std::cout << "  Beamspot size: " << beamspot << " m\n";
	std::cout << "  Beam focal point: " << focalpoint << " m\n";
	std::cout << "  Total size on target: " << detectorspot << " m\n";
	
	// Setup the particle range tables
	for(unsigned int i = 0; i < 6; i++){
		std::cout << " Setting up range table for '" << particles[i].GetName() << "' in " << BC408.GetName() << "...";
		if(particles[i].SetMaterial(&BC408, Ebeam0, EbeamSpread)){ std::cout << " done\n"; }
		else{ std::cout << " failed\n"; }
	}

	// Root stuff
	TFile *file = new TFile("PHOSWICH.root", "RECREATE");
	TTree *VIKARtree = new TTree("VIKAR", "VIKAR output tree");
	
	RecoilObject RECOILdata;
	VIKARtree->Branch("Recoil", &RECOILdata);
	
	TH2D *hist = new TH2D("hist", "Vandle Phoswich dE vs. E", 400, 0, 35, 400, 0, 35);
	TCanvas *can = new TCanvas("can", "Phoswich");
	can->cd();
	
	hist->SetStats(false);
	hist->GetXaxis()->SetTitle("E (MeV)");
	hist->GetYaxis()->SetTitle("dE (MeV)");
	
	// Do a monte carlo simulation
	bool hit;
	int face1, face2;
	double Ebeam, hit_x, hit_y, hit_z;
	Vector3 lab_beam_spot, lab_beam_direction;
	Vector3 lab_beam_start(0.0, 0.0, focalpoint);
	Vector3 HitDetect1, HitDetect2;
	Vector3 path, stop_point;
	Vector3 detect_sphere;
	
	double range, path_length;
	double tof, Enew, qdc1, qdc2;
	unsigned int count = 0;
	unsigned int total_count = 0;
	for(unsigned int i = 0; i < 10000; i++){
		for(unsigned int j = 0; j < 6; j++){
			qdc1 = 0.0; qdc2 = 0.0;
			RandomCircle(beamspot, 0.0, lab_beam_spot); 
			lab_beam_direction = lab_beam_spot-lab_beam_start;
		
			count = 0;
			Ebeam = Ebeam0 + rndgauss0(EbeamSpread); 
			for(std::vector<Planar*>::iterator iter = detectors.begin(); iter != detectors.end(); iter++){
				hit = (*iter)->IntersectPrimitive(lab_beam_start, lab_beam_direction, HitDetect1, HitDetect2, face1, face2, hit_x, hit_y, hit_z);
				if(hit){ 
					path = HitDetect2-HitDetect1;
					path_length = path.Length();
					range = particles[j].GetRange(Ebeam);
					if(range < path_length){ // The particle stops in the detector (E)
						stop_point = HitDetect1 + path*range;
						Cart2Sphere(stop_point, detect_sphere);
						qdc1 = Ebeam;
						tof = range*std::sqrt(particles[j].GetA()/(2*qdc1*1.60217657E-13*6.02214129E26)) + rndgauss0(timeRes);
						RECOILdata.Append(stop_point.axis[0], stop_point.axis[1], stop_point.axis[2], detect_sphere.axis[1]*rad2deg,
										 detect_sphere.axis[2]*rad2deg, qdc1, tof*(1E9), hit_x, hit_y, hit_z, count);
						break;
					}
					else{ // The particle passes through, but loses some energy (dE)
						Cart2Sphere(HitDetect2, detect_sphere);
						qdc2 = particles[j].GetEnergy(path_length);
						tof = path_length*std::sqrt(particles[j].GetA()/(2*qdc2*1.60217657E-13*6.02214129E26)) + rndgauss0(timeRes);
						RECOILdata.Append(HitDetect2.axis[0], HitDetect2.axis[1], HitDetect2.axis[2], detect_sphere.axis[1]*rad2deg,
										 detect_sphere.axis[2]*rad2deg, qdc2, tof*(1E9), hit_x, hit_y, hit_z, count);
						Ebeam = Ebeam - qdc2;
						total_count++;
					}
				}
				count++;
			}
			hist->Fill(qdc1, qdc2);
		}
		if(RECOILdata.recoil_mult > 0){ VIKARtree->Fill(); }
		RECOILdata.Zero();
	}
	
	VIKARtree->Write();
	std::cout << " Done! Detected " << total_count << " hits (" << 100.0*total_count/60000.0 << "%)\n";
	std::cout << "  Wrote file 'PHOSWICH.root'\n";
	std::cout << "   Wrote " << VIKARtree->GetEntries() << " tree entries for VIKAR\n";

	hist->Draw("COLZ");
	can->WaitPrimitive();

	file->Close();
	delete file;	

	return 0;
}
