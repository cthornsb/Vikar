#include <iostream>
#include <fstream>
#include <vector>

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

#define ANTHRACENE 17400 // Anthracene photon output (1/MeV)

// Fast pulse parameters
const double beta_fast = 1.7750575; // Decay parameter of the pulse exponential in ns
const double gamma_fast = 115.64125; // Width of the inverted square gaussian of the pulse in ns^4
const double coeff_fast = 0.371821719; // Normalization parameter
const double LO_fast = 0.65*ANTHRACENE; // Light output (1/MeV)

// Slow pulse parameters
const double beta_slow = 125.0; // Decay parameter of the pulse exponential in ns
const double gamma_slow = 3.0E4; // Width of the inverted square gaussian of the pulse in ns^4
const double coeff_slow = 112.6387378; // Normalization parameter
const double LO_slow = 0.41*ANTHRACENE; // Light output (1/MeV)

const double pixie_time_res = 4.0; // Pixie-16 time resolution (4 ns for 250 MHz) in ns
const double pixie_pulse_width = 1000.0; // Width of Pixie-16 pulse window in ns
const unsigned int num_samples = (unsigned int)(pixie_pulse_width/pixie_time_res);

/*struct BirksTable{
	unsigned int num_entries;
	double *energy, *light;
	double L0, kB, C;
	bool use_table;
	
	BirksTable(){
		num_entries = 0; L0 = 0.0; kB = 0.0; C = 0.0;
		energy = NULL; light = NULL;
		use_table = false;
	}
	
	bool Init(unsigned int num_entries_, double L0_, double kB_, double C_, double startE_, double stopE_, double Z_, double mass_, Material *mat_){
		if(use_table){ return false; }
		L0 = L0_;
		kB = kB_;
		C = C_;
		
		// Use Material to fill the arrays
		num_entries = num_entries_;
		double step = (stopE_-startE_)/(num_entries_-1);
		for(unsigned int i = 0; i < num_entries_; i++){
			energy[i] = startE_ + i*step; // Energy in MeV
			light[i] = mat_->Birks(energy[i], Z_, mass_); // Stopping power in MeV/m
		}
		use_table = true;

		return true;
	}
	
	double GetLightOutput(double energy_){
		if(!use_table){ return -1; }
		for(unsigned int i = 0; i < num_entries-1; i++){
			if(energy_ == energy[i]){ return range[i]; }
			else if(energy_ == energy[i+1]){ return range[i+1]; }
			else if(energy_ >= energy[i] && energy_ <= energy[i+1]){
				// Interpolate and return the result
				return (((light[i+1]-light[i])/(energy[i+1]-energy[i]))*(energy_-energy[i])+light[i]);
			}
		}
		return -1;
	}
};*/

// Fast pulse fitting function
//  t_: Time in ns
//  alpha_: normalization of pulse
//  phi_: phase of pulse 1 in ns
double fast_pulse_function(double t_, double alpha_, double phi_){
	double arg = t_ - phi_;
	if(arg >= 0.0){ return alpha_*std::exp(-arg/beta_fast)*(1 - std::exp(-arg*arg*arg*arg/gamma_fast)); }
	return 0.0;
}

// Slow pulse fitting function
//  t_: Time in ns
//  alpha_: normalization of pulse
//  phi_: phase of pulse 1 in ns
double slow_pulse_function(double t_, double alpha_, double phi_){
	double arg = t_ - phi_;
	if(arg >= 0.0){ return alpha_*std::exp(-arg/beta_slow)*(1 - std::exp(-arg*arg*arg*arg/gamma_slow)); }
	return 0.0;
}

void discretize(double Afast_, double Aslow_, double phi_, std::vector<int> &pulse){
	pulse.clear();
	double fast_component, slow_component;
	for(unsigned int i = 0; i <= num_samples; i++){
		fast_component = fast_pulse_function(i*pixie_time_res, Afast_, phi_);
		slow_component = slow_pulse_function(i*pixie_time_res, Aslow_, phi_);
		pulse.push_back((int)(fast_component + slow_component));
	}
}

double calc_mass(double Z_, double A_, double BE_A_=0.0){
	return Z_*proton_RME + (A_-Z_)*neutron_RME - A_*BE_A_;
}

int main(int argc, char *argv[]){
	bool use_focus = true;
	bool use_pulse = false;
	if(argc > 1){
		if(strcmp(argv[1], "debug") == 0){
			std::cout << " Debugging! Using cylindrical beam...\n";
			use_focus = false;
		}
		else if(strcmp(argv[1], "pulse") == 0){
			std::cout << " Saving phoswich pulses to file...\n";
			use_pulse = true;
		}
	}

	char* dummy[0]; 
	gSystem->Load("libTree");
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	
	double Ebeam0 = 27.0; // MeV
	double EbeamSpread = 1.6; // MeV
	double thetaMax = 13.3241*deg2rad; // radians
	double beamspot = 0.02; // m
	
	double detectorspot = beamspot;
	double focalpoint = 0.0;
	if(use_focus){
		detectorspot = beamspot+0.3*std::tan(thetaMax); // m
		focalpoint = -beamspot/(2.0*std::tan(thetaMax)); // m
	}
	
	std::vector<Planar*> detectors;
	unsigned int num_dets = ReadDetFile("./detectors/Phoswich.det", detectors);
	if(num_dets == 0){
		std::cout << " Error: Failed to load detector file\n";
		return 1;
	}
	
	std::cout << " Successfully loaded " << num_dets << " detectors from file\n";

	// All of the possible beam particles
	Particle particles[6];
	particles[0].SetParticle("Li6", 3, 6, 5.332331); // Lithium-6
	particles[1].SetParticle("Li7", 3, 7, 5.606439); // Lithium-7
	particles[2].SetParticle("Be7", 4, 7, 5.371548); // Beryillium-7
	particles[3].SetParticle("Be8", 4, 8, 7.062435); // Beryillium-8 (short lived)
	particles[4].SetParticle("B8", 5, 8, 7.062435); // Boron-8
	particles[5].SetParticle("B9", 5, 9, 6.257070); // Boron-9 (short lived)
	
	Material BC408;
	BC408.ReadMatFile("./materials/BC408.mat");
	std::cout << "  Beam energy: " << Ebeam0 << " MeV\n";
	std::cout << "  Energy spread: " << EbeamSpread << " MeV\n";
	std::cout << "  Theta max: " << thetaMax << " rad\n";
	std::cout << "  Beam size on target: " << beamspot << " m\n";
	if(use_focus){ 
		std::cout << "  Beam focal point: " << focalpoint << " m\n"; 
		std::cout << "  Total size on detector: " << detectorspot << " m\n";
	}
	
	//BirksTable tables[6];
	//tables[0].Init(100, LO_fast, ; // Lithium-6 in BC408
	
	// Setup the particle range tables
	for(unsigned int i = 0; i < 6; i++){
		std::cout << " Setting up range table for '" << particles[i].GetName() << "' in BC408...";
		if(particles[i].SetMaterial(&BC408, Ebeam0, EbeamSpread)){ std::cout << " done\n"; }
		else{ std::cout << " failed\n"; }
	}

	// Setup the slow scintillator light output tables
	/*for(unsigned int i = 0; i < 6; i++){
		std::cout << " Setting up light output table for '" << particles[i].GetName() << "' in BC408...";
		if(tables[i].Init(100, LO_slow, )){ std::cout << " done\n"; }
		else{ std::cout << " failed\n"; }
	}*/

	// Root stuff
	TFile *file = new TFile("PHOSWICH.root", "RECREATE");
	TTree *VIKARtree = new TTree("VIKAR", "VIKAR output tree");
	
	RecoilObject RECOILdata;
	PulseObject PULSEdata;
	std::vector<int> pulse;
	
	if(!use_pulse){ VIKARtree->Branch("Recoil", &RECOILdata); }
	else{ VIKARtree->Branch("Pulse", &PULSEdata); }
	
	TH2D *hist = new TH2D("hist", "Vandle Phoswich dE vs. E", 100, 0, 35, 100, 0, 35);
	TCanvas *can = new TCanvas("can", "Phoswich");
	can->cd();
	
	hist->SetStats(false);
	hist->GetXaxis()->SetTitle("E (MeV)");
	hist->GetYaxis()->SetTitle("dE (MeV)");
	
	// Do a monte carlo simulation
	bool hit;
	int face1, face2;
	double Ebeam, hit_x, hit_y, hit_z;
	Vector3 lab_beam_spot;

	Vector3 HitDetect1, HitDetect2;
	Vector3 path, stop_point;
	Vector3 detect_sphere;
	
	Vector3 lab_beam_start(0.0, 0.0, focalpoint);	
	Vector3 lab_beam_direction(0.0, 0.0, 1.0);
	
	double range, path_length;
	double qdc1, qdc2;
	unsigned int count = 0;
	unsigned int total_count = 0;
	for(unsigned int i = 0; i < 100; i++){
		for(unsigned int j = 0; j < 6; j++){
			qdc1 = 0.0; qdc2 = 0.0;
			RandomCircleUp(beamspot, 0.0, lab_beam_spot); 
			if(use_focus){ lab_beam_direction = lab_beam_spot-lab_beam_start; }
			else{ lab_beam_start = lab_beam_spot - Vector3(0.0, 0.0, 1.0); }
		
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
						
						if(!use_pulse){ 
							RECOILdata.Append(stop_point.axis[0], stop_point.axis[1], stop_point.axis[2], detect_sphere.axis[1]*rad2deg,
											  detect_sphere.axis[2]*rad2deg, qdc1, j, hit_x, hit_y, hit_z, count);
						}
												
						break;
					}
					else{ // The particle passes through, but loses some energy (dE)
						Cart2Sphere(HitDetect2, detect_sphere);
						qdc2 = Ebeam - particles[j].GetNewE(Ebeam, path_length);
						
						if(!use_pulse){ 
							RECOILdata.Append(HitDetect2.axis[0], HitDetect2.axis[1], HitDetect2.axis[2], detect_sphere.axis[1]*rad2deg,
											 detect_sphere.axis[2]*rad2deg, qdc2, j, hit_x, hit_y, hit_z, count);
						}
						
						Ebeam = Ebeam - qdc2;
						total_count++;
					}
				}
				count++;
			}
			if(qdc1 != 0.0 && qdc2 != 0.0){ 
				hist->Fill(qdc1, qdc2); 
				if(use_pulse){
					// Use the fast/slow scintillator properties to construct pulses from the energies
					//discretize(LO_fast*qdc2/coeff_fast, LO_slow*qdc1/coeff_slow, 100.0, pulse);
					discretize(0.0, LO_slow*qdc1/coeff_slow, 100.0, pulse);
					//discretize(LO_fast*qdc2/coeff_fast, 0.0, 100.0, pulse);
					
					// Copy the vector and fill the tree
					PULSEdata.Set(pulse, j);
					VIKARtree->Fill();
					PULSEdata.Zero();
				}
			}
		}
		if(!use_pulse){
			if(RECOILdata.recoil_mult > 0){ VIKARtree->Fill(); }
			RECOILdata.Zero();
		}

	}
	
	VIKARtree->Write();
	std::cout << " Done! Detected " << total_count << " hits (" << 100.0*total_count/60000.0 << "%)\n";
	std::cout << "  Wrote file 'PHOSWICH.root'\n";
	std::cout << "   Wrote " << VIKARtree->GetEntries() << " tree entries for VIKAR\n";

	hist->Draw("COLZ");
	can->WaitPrimitive();

	file->Close();
	delete file;	

	rootapp->Delete();

	return 0;
}
