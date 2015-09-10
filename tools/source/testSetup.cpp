#include <iostream>

#include "vikar_core.h"
#include "detectors.h"
#include "Structures.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

struct DataPack{
	TFile *file;
	TTree *tree;
	TH1D *hist;
	bool init;

	MonteCarloStructure MCARLOdata;
	ReactionObjectStructure REACTIONdata;

	DataPack(){
		file = NULL;
		tree = NULL;
		hist = NULL;
		init = false;
	}
	
	DataPack(std::string fname_, unsigned int Nbins_, double low_, double high_, bool write_rxn_){
		Open(fname_, Nbins_, low_, high_, write_rxn_);
	}

	~DataPack(){
		if(init){ Close(); }
	}

	bool IsInit(){ return init; }

	bool Open(std::string fname_, unsigned int Nbins_, double low_, double high_, bool write_rxn_){
		if(init){ return false; }

		file = new TFile(fname_.c_str(), "RECREATE");
		
		if(!file->IsOpen()){
			init = false;
			delete file;
		}
		
		tree = new TTree("VIKAR", "Monte carlo detector efficiency tree");
		hist = new TH1D("hist", "Detector Hits vs. CoM Angle", Nbins_, low_, high_);
		
		tree->Branch("MCarlo", &MCARLOdata);
		if(write_rxn_){ tree->Branch("Reaction", &REACTIONdata); }
		
		return (init = true);
	}

	bool Close(){
		if(!init){ return false; }
		
		file->cd();
		tree->Write();
		hist->Write();
		file->Close();
		
		init = false;
		
		std::cout << "  Wrote monte carlo file 'mcarlo.root'\n";
		
		return true;
	}
	
	void Zero(){
		MCARLOdata.Zero();
		REACTIONdata.Zero();
	}
};

double proper_value(const std::string &prompt_, const double &min_=0.0, bool ge_=false){
	double output = -1;
	while(true){
		std::cout << " " << prompt_;
		std::cin >> output;
		if(!ge_){
			if(output > min_){ break; }
			std::cout << "  Error: Invalid value! Input must be > " << min_ << ".\n";
		}
		else{
			if(output >= min_){ break; }
			std::cout << "  Error: Invalid value! Input must be >= " << min_ << ".\n";
		}
	}
	return output;
}

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates one output root file named 'mcarlo.root'
// fwhm_ (m) allows the use of a gaussian particle "source". If fwhm_ == 0.0, a point source is used
// angle_ (rad) allows the rotation of the particle source about the y-axis
unsigned int TestDetSetup(DataPack *pack, Planar *bar_array, unsigned int num_bars, unsigned int num_trials, bool WriteRXN_, double fwhm_, double angle_, bool ejectile_=true){	
	if(!pack || !bar_array){ return 0; }
	double tempx, tempy, tempz;
	double dummyR, hitTheta, hitPhi;
	unsigned int count, total;
	unsigned int bar;
	Vector3 flight_path;
	Vector3 temp_vector1;
	Vector3 temp_vector2;
	Vector3 temp_ray, offset, dummyVector;
	int face1, face2, type;
	bool found_hit;
	
	Matrix3 matrix;
	bool use_gaussian_beam = true;
	bool use_rotated_source = false;
	if(fwhm_ == 0.0){ 
		offset = zero_vector; 
		use_gaussian_beam = false;
	}
	if(angle_ != 0.0){ 
		use_rotated_source = true; 
		matrix.SetRotationMatrixSphere(angle_, 0.0);
	}
	
	total = 0; count = 0; type = 0;
	while(count < num_trials){
		found_hit = false;
	
		UnitSphereRandom(temp_ray); // Generate a uniformly distributed random point on the unit sphere
		if(use_gaussian_beam){ // Generate an offset based on a gaussian beam
			double x_offset = rndgauss0(fwhm_);
			double y_offset = rndgauss0(fwhm_);
			offset.axis[0] = x_offset;
			offset.axis[1] = y_offset;
			offset.axis[2] = 0.0;
			if(use_rotated_source){ matrix.Transform(offset); } // This will rotate the "source" about the y-axis
		}
		if(WriteRXN_){ 
			pack->REACTIONdata.Zero();
			pack->REACTIONdata.Append(0.0, 0.0, offset.axis[0], offset.axis[1], offset.axis[2], temp_ray.axis[0], temp_ray.axis[1], temp_ray.axis[2]); 
		}
		for(bar = 0; bar < num_bars; bar++){
			if(bar_array[bar].IntersectPrimitive(offset, temp_ray, temp_vector1, temp_vector2, face1, face2, tempx, tempy, tempz)){
				if(bar_array[bar].IsEjectileDet()){
					if(bar_array[bar].IsRecoilDet()){ type = 2; } // Both ejectile & recoil
					else{ type = 0; } // Ejectile
				}
				else if(bar_array[bar].IsRecoilDet()){ type = 1; } // Recoil
				else{ continue; }
				
				dummyVector = (offset + temp_ray);
				Cart2Sphere(dummyVector.axis[0], dummyVector.axis[1], dummyVector.axis[2], dummyR, hitTheta, hitPhi);
				
				// Is this a valid event?
				if(type == 2 || (ejectile_ && type == 0) || (!ejectile_ && type == 1)){
					pack->MCARLOdata.Append(temp_vector1.axis[0], temp_vector1.axis[1], temp_vector1.axis[2],
										   temp_vector2.axis[0], temp_vector2.axis[1], temp_vector2.axis[2],
										   hitTheta*rad2deg, hitPhi*rad2deg, face1, face2, bar, type);
					found_hit = true;
				}
			}
		}
		
		if(WriteRXN_ || found_hit){
			pack->tree->Fill();
			pack->REACTIONdata.Zero();
			pack->MCARLOdata.Zero();
			pack->MCARLOdata.Zero();
			if(found_hit){ count++; }
		}
		
		total++;
	}

	return total;
}

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates an output ascii file named 'mcarlo.hist' which contains 
// bin contents of a histogram of event center of mass angles.
// fwhm_ (m) allows the use of a gaussian particle "source". If fwhm_ == 0.0, a point source is used
// angle_ (rad) allows the rotation of the particle source about the y-axis
unsigned int TestDetSetup(DataPack *pack, Planar *bar_array, Kindeux *kind_, unsigned int num_bars, unsigned int num_trials, double beamE_,
						  unsigned int num_bins_, double start_, double stop_, double fwhm_/*=0.0*/, double angle_/*=0.0*/){
	if(!pack || !bar_array || !kind_ || !kind_->IsInit()){ return 0; }
	double tempx, tempy, tempz;
	unsigned int count, total;
	unsigned int bar;
	Vector3 flight_path;
	Vector3 temp_vector1;
	Vector3 temp_vector2;
	Vector3 Ejectile, EjectSphere;
	Vector3 Recoil, RecoilSphere;
	Vector3 offset, dummyVector;
	int face1, face2;
	
	double ejectE, recoilE;
	double com_angle;
	bool eject_hit, recoil_hit;
	
	Matrix3 matrix;
	bool use_gaussian_beam = true;
	bool use_rotated_source = false;
	if(fwhm_ == 0.0){ 
		offset = zero_vector; 
		use_gaussian_beam = false;
	}
	if(angle_ != 0.0){ 
		use_rotated_source = true; 
		matrix.SetRotationMatrixSphere(angle_, 0.0);
	}

	double dummyR, dummyPhi, ejectTheta, recoilTheta, alpha;
	double m1m2 = kind_->GetMbeam()/kind_->GetMtarg();

	total = 0; count = 0;
	while(count < num_trials){
		if(use_gaussian_beam){ // Generate an offset based on a gaussian beam
			double x_offset = rndgauss0(fwhm_);
			double y_offset = rndgauss0(fwhm_);
			offset.axis[0] = x_offset;
			offset.axis[1] = y_offset;
			offset.axis[2] = 0.0;
			if(use_rotated_source){ matrix.Transform(offset); } // This will rotate the "source" about the y-axis
		}
		
		kind_->FillVars(beamE_, ejectE, recoilE, EjectSphere, RecoilSphere, com_angle);
		Sphere2Cart(EjectSphere, Ejectile); 
		Sphere2Cart(RecoilSphere, Recoil);
		
		eject_hit = false; recoil_hit = false;
		for(bar = 0; bar < num_bars; bar++){
			if(!eject_hit && bar_array[bar].IsEjectileDet()){
				if(bar_array[bar].IntersectPrimitive(offset, Ejectile, temp_vector1, dummyVector, face1, face2, tempx, tempy, tempz)){ eject_hit = true; }
			}
			else if(!recoil_hit && bar_array[bar].IsRecoilDet()){
				if(bar_array[bar].IntersectPrimitive(offset, Recoil, temp_vector2, dummyVector, face1, face2, tempx, tempy, tempz)){ recoil_hit = true; }
			}
			
			if(eject_hit && recoil_hit){
				// Encountered coincidence event
				//com_angle = Ejectile.Dot(Recoil)*rad2deg;
				
				Cart2Sphere(Ejectile, dummyR, ejectTheta, dummyPhi);
				Cart2Sphere(Recoil, dummyR, recoilTheta, dummyPhi);
				
				alpha = std::tan(recoilTheta)/std::tan(ejectTheta);
				com_angle = std::acos((1-alpha*m1m2)/(1+alpha))*rad2deg;
				
				pack->hist->Fill(com_angle);

				count++;
				break;
			}
		}
		total++;
	}
	
	return total;
}

void help(char * prog_name_){
	std::cout << "  SYNTAX: " << prog_name_ << " [detfile]\n";
}

int main(int argc, char *argv[]){
	if(argc < 2){
		std::cout << " Error: Invalid number of arguments to " << argv[0] << ". Expected 1, received " << argc-1 << ".\n";
		help(argv[0]);
		return 1;
	}

	bool do_hist_run = false;

	Planar *detectors = NULL;

	std::cout << " Reading in NewVIKAR detector setup file...\n";
	int Ndet = ReadDetFile(argv[1], detectors);
	if(Ndet < 0){
		std::cout << " Error: failed to load detector setup file!\n";
		return 1; 
	}
	else if(Ndet == 0){ std::cout << " Error: Found no detectors in the detector setup file!\n"; }

	unsigned int Nrecoil = 0;
	unsigned int Nejectile = 0;
	for(int index = 0; index < Ndet; index++){
		if(detectors[index].IsRecoilDet()){ Nrecoil++; }
		if(detectors[index].IsEjectileDet()){ Nejectile++; }
	}

	// Check there's at least 1 detector!
	if((Nejectile+Nrecoil) < 1){
		std::cout << " Error: Found no valid detectors in the detector setup file!\n"; 
		return 1;
	}

	// Report on how many detectors were read in
	std::cout << " Found " << (Nejectile+Nrecoil) << " detectors in file " << argv[1] << "\n";
	std::cout << "  " << Nejectile << " ejectile detectors and " << Nrecoil << " recoil detectors\n\n";

	double beamspot = 0.0, targangle = 0.0;
	unsigned int Nwanted = 0;
	unsigned int total_found = 0;
	bool WriteReaction = false;

	beamspot = proper_value("Enter source FWHM (m) (0 for point source): ", 0.0, true);
	if(beamspot > 0.0){ targangle = proper_value("Enter target angle (degrees): ", 0.0, true); }
	std::cout << " Write reaction data? "; std::cin >> WriteReaction; 

	DataPack pack;

	unsigned int Nbins;
	double startAngle;
	double stopAngle;

	if(do_hist_run){
		Nbins = (unsigned int)proper_value("Enter number of CoM bins: ", 0.0, false);
		startAngle = proper_value("Enter start CoM angle (degrees): ", 0.0, true);
		stopAngle = proper_value("Enter stop CoM angle (degrees): ", startAngle, false);
		pack.Open("mcarlo.root", Nbins, startAngle, stopAngle, WriteReaction);
	}
	else{ pack.Open("mcarlo.root", 1, 0, 1, WriteReaction); }

	std::cout << std::endl;

	if(Nejectile > 0){		
		Nwanted = (unsigned int)proper_value("Enter number of ejectile MC events: ", 0.0, true);

		// Process ejectile detectors
		if(Nwanted > 0){
			std::cout << "  Performing Monte Carlo test on ejectile detectors...\n";
			total_found = TestDetSetup(&pack, detectors, Ndet, Nwanted, WriteReaction, beamspot, targangle, true);
			if(do_hist_run){
				//total_found = TestDetSetup(&pack, detectors, &kind, Ndet, Nwanted, Ebeam, Nbins, startAngle*deg2rad, stopAngle*deg2rad, beamspot, targangle));
			}
		
			std::cout << "  Found " << Nwanted << " ejectile events in " << total_found << " trials (" << 100.0*Nwanted/total_found << "%)\n\n";
		}
	}
	else{ std::cout << " Note: Found no ejectile detectors in detector file.\n\n"; }

	if(Nrecoil > 0){
		Nwanted = (unsigned int)proper_value("Enter number of recoil MC events: ", 0.0, true);

		// Process recoil detectors
		if(Nwanted > 0){
			std::cout << "  Performing Monte Carlo test on recoil detectors...\n";
			total_found = TestDetSetup(&pack, detectors, Ndet, Nwanted, WriteReaction, beamspot, targangle, false);
			if(do_hist_run){
				//total_found = TestDetSetup(&pack, detectors, &kind, Ndet, Nwanted, Ebeam, Nbins, startAngle*deg2rad, stopAngle*deg2rad, beamspot, targangle));
			}
		
			std::cout << "  Found " << Nwanted << " recoil events in " << total_found << " trials (" << 100.0*Nwanted/total_found << "%)\n\n";
		}
	}
	else{ std::cout << " Note: Found no recoil detectors in detector file.\n\n"; }

	std::cout << " Finished geometric efficiency test on detector setup...\n";
	
	delete[] detectors;

	return 0;
}
