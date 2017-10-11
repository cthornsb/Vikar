#include <iostream>
#include <sstream>
#include <fstream>

#include "vandmc_core.h"
#include "kindeux.h"
#include "detectors.h"
#include "materials.h"
#include "Structures.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TNamed.h"

Vector3 zero_vector(0.0, 0.0, 0.0);

struct DataPack{
	TFile *file;
	TTree *tree;
	bool init;

	double offsetX, offsetY, offsetZ;
	double trajX, trajY, trajZ;
	double Ereact, Eeject, Erecoil;
	double comAngle;

	MonteCarloStructure MCARLOdata;

	reactData reaction;

	DataPack(){
		file = NULL;
		tree = NULL;
		init = false;
	}
	
	DataPack(std::string fname_, unsigned int Nbins_, bool write_rxn_=false){
		Open(fname_, write_rxn_);
	}

	~DataPack(){
		if(init){ Close(); }
	}

	bool IsInit(){ return init; }

	bool Open(std::string fname_, bool write_rxn_){
		if(init){ return false; }

		file = new TFile(fname_.c_str(), "RECREATE");
		
		if(!file->IsOpen()){
			init = false;
			delete file;
		}
		
		tree = new TTree("data", "Monte carlo detector efficiency tree");
		
		tree->Branch("mcarlo", &MCARLOdata);
		if(write_rxn_){ 
			tree->Branch("Ereact", &Ereact); 
			tree->Branch("Eeject", &Eeject); 
			tree->Branch("Erecoil", &Erecoil); 
			tree->Branch("comAngle", &comAngle); 
		}
		
		return (init = true);
	}

	bool Close(){
		if(!init){ return false; }
		
		file->cd();
		tree->Write();
		file->Close();
		
		init = false;
		
		std::cout << "  Wrote monte carlo file 'mcarlo.root'\n";
		
		return true;
	}
	
	void Zero(){
		MCARLOdata.Zero();
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
unsigned int TestDetSetup(DataPack *pack, const std::vector<Primitive*> &bar_array, unsigned int num_trials, bool WriteRXN_, double fwhm_, double angle_, bool ejectile_=true, Kindeux *kind=NULL){
	if(!pack){ return 0; }
	double dummyR, hitTheta, hitPhi;
	double t1, t2;
	unsigned int count, total;
	Vector3 flight_path;
	Vector3 temp_vector;
	Vector3 temp_normal;
	Vector3 temp_ray, offset, dummyVector;
	Vector3 ejectileVector, recoilVector;
	int type;
	bool found_hit;
	
	bool useKindeux = (kind != 0x0);
	
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
	
	unsigned int num_trials_chunk = num_trials/10;
	unsigned int chunk_num = 1;
	
	total = 0; count = 0; type = 0;
	while(count < num_trials){
		if(count != 0 && count == num_trials_chunk*chunk_num){ // Print a status update.
			std::cout << "  " << (chunk_num++)*10 << "% - " << "Detected " << count << " of " << total << " total events (" << count*100.0/total << "%)\n";
		}
	
		found_hit = false;

		if(useKindeux){ // Generate a uniformly distributed random point on the unit sphere in the center-of-mass frame.
			if(ejectile_){
				kind->FillVars(pack->reaction, ejectileVector, dummyVector);
				Sphere2Cart(ejectileVector, temp_ray);
			}
			else{
				kind->FillVars(pack->reaction, dummyVector, recoilVector);
				Sphere2Cart(recoilVector, temp_ray);
			}
		}
		else{ // Generate a uniformly distributed random point on the unit sphere in the lab frame.
			UnitSphereRandom(temp_ray);
		}

		// Generate an offset based on a gaussian beam.
		if(use_gaussian_beam){
			double x_offset = rndgauss0(fwhm_);
			double y_offset = rndgauss0(fwhm_);
			offset.axis[0] = x_offset;
			offset.axis[1] = y_offset;
			offset.axis[2] = 0.0;
			if(use_rotated_source){ matrix.Transform(offset); } // This will rotate the "source" about the y-axis
		}

		for(std::vector<Primitive*>::const_iterator iter = bar_array.begin(); iter != bar_array.end(); iter++){
			if((*iter)->IsEjectileDet()){
				if((*iter)->IsRecoilDet()){ type = 2; } // Both ejectile & recoil
				else{ type = 0; } // Ejectile
			}
			else if((*iter)->IsRecoilDet()){ type = 1; } // Recoil
			else{ continue; }
		
			// Is this detector valid for this type of event?
			if(ejectile_ && type == 1){ continue; }
			else if(!ejectile_ && type == 0){ continue; }
			
			// Check for intersections with detectors.
			if((*iter)->IntersectPrimitive(offset, temp_ray, temp_vector, temp_normal, t1, t2)){
				if(t2 > 0){ dummyVector = (offset + temp_ray * t2); }
				else{ dummyVector = Vector3(0,0,0); }
				
				Cart2Sphere(temp_vector, dummyR, hitTheta, hitPhi);
				
				pack->MCARLOdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2],
									   dummyVector.axis[0], dummyVector.axis[1], dummyVector.axis[2],
									   hitTheta*rad2deg, hitPhi*rad2deg, 0, 0, (*iter)->GetLoc(), type);
									   
				found_hit = true;
			}
		}
		
		if(WriteRXN_){
			pack->Ereact = pack->reaction.Ereact;
			pack->Eeject = pack->reaction.Eeject;
			pack->Erecoil = pack->reaction.Erecoil;
			pack->comAngle = pack->reaction.comAngle*rad2deg;
		}

		if(found_hit){
			pack->tree->Fill();
			pack->MCARLOdata.Zero();
			pack->MCARLOdata.Zero();
			if(found_hit){ count++; }
		}
		else if(WriteRXN_){
			pack->tree->Fill();
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

	std::vector<Primitive*> detectors;

	std::cout << " Reading in NewVIKAR detector setup file...\n";
	int Ndet = ReadDetFile(argv[1], detectors);
	if(Ndet < 0){
		std::cout << " Error: failed to load detector setup file!\n";
		return 1; 
	}
	else if(Ndet == 0){ std::cout << " Error: Found no detectors in the detector setup file!\n"; }

	unsigned int Nrecoil = 0;
	unsigned int Nejectile = 0;
	for(std::vector<Primitive*>::iterator iter = detectors.begin(); iter != detectors.end(); iter++){
		if((*iter)->IsRecoilDet()){ Nrecoil++; }
		if((*iter)->IsEjectileDet()){ Nejectile++; }
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
	bool UseKinematics = false;

	beamspot = proper_value("Enter source FWHM (m) (0 for point source): ", 0.0, true);
	if(beamspot > 0.0){ targangle = proper_value("Enter target angle (degrees): ", 0.0, true); }
	std::cout << " Write reaction data? "; std::cin >> WriteReaction; 
	std::cout << " Use kinematics? "; std::cin >> UseKinematics;
	
	DataPack pack;
	
	Kindeux *kind = NULL;
	if(UseKinematics){
		double Mbeam[3], Mtarg[3];
		double Mrecoil[4], Meject[4];
		double recoilStates[1];
		double Qgs;
		bool userChoice;

		Particle beam, targ, recoil, ejectile;

		std::cout << "  Load reaction parameters from file? "; std::cin >> userChoice;

		std::string rxnFilename;
		if(userChoice){
			std::cout << "  Enter input filename: "; std::cin >> rxnFilename;
			std::ifstream rxnFile(rxnFilename.c_str());
			rxnFile >> Mbeam[0] >> Mbeam[1] >> Mbeam[2];
			rxnFile >> Mtarg[0] >> Mtarg[1] >> Mtarg[2];
			rxnFile >> Mrecoil[0] >> Mrecoil[1] >> Mrecoil[2] >> Mrecoil[3];
			rxnFile >> Meject[0] >> Meject[1] >> Meject[2] >> Meject[3];
			rxnFile >> pack.reaction.Ereact;
			rxnFile.close();
		}
		else{
			std::cout << "  Enter beam Z, A, BE/A(MeV): "; std::cin >> Mbeam[0] >> Mbeam[1] >> Mbeam[2];
			std::cout << "  Enter target Z, A, BE/A(MeV): "; std::cin >> Mtarg[0] >> Mtarg[1] >> Mtarg[2];
			std::cout << "  Enter recoil Z, A, BE/A(MeV), Excitation(MeV): "; std::cin >> Mrecoil[0] >> Mrecoil[1] >> Mrecoil[2] >> Mrecoil[3];
			std::cout << "  Enter ejectile Z, A, BE/A(MeV), Excitation(MeV): "; std::cin >> Meject[0] >> Meject[1] >> Meject[2] >> Meject[3];
			std::cout << "  Enter beam energy (MeV): "; std::cin >> pack.reaction.Ereact;
			std::cout << "  Write reaction parameters to file? "; std::cin >> userChoice;
			if(userChoice){
				std::cout << "  Enter output filename: "; std::cin >> rxnFilename;
				std::ofstream rxnFile(rxnFilename.c_str());
				rxnFile << Mbeam[0] << "\t" << Mbeam[1] << "\t" << Mbeam[2] << "\n";
				rxnFile << Mtarg[0] << "\t" << Mtarg[1] << "\t" << Mtarg[2] << "\n";
				rxnFile << Mrecoil[0] << "\t" << Mrecoil[1] << "\t" << Mrecoil[2] << "\t" << Mrecoil[3] << "\n";
				rxnFile << Meject[0] << "\t" << Meject[1] << "\t" << Meject[2] << "\t" << Meject[3] << "\n";
				rxnFile << pack.reaction.Ereact;
				rxnFile.close();
			}
		}

		beam.SetParticle("", Mbeam[0], Mbeam[1], Mbeam[2]);
		targ.SetParticle("", Mtarg[0], Mtarg[1], Mtarg[2]);
		recoil.SetParticle("", Mrecoil[0], Mrecoil[1], Mrecoil[2]);
		ejectile.SetParticle("", Meject[0], Meject[1], Meject[2]);

		recoilStates[0] = Mrecoil[3];

		Qgs = (beam.GetMass()+targ.GetMass())-(recoil.GetMass()+ejectile.GetMass());

		std::cout << "\n Beam energy (MeV)       = " << pack.reaction.Ereact << std::endl;
		std::cout << " Beam mass (MeV/c^2)     = " << beam.GetMass() << std::endl;
		std::cout << " Target mass (MeV/c^2)   = " << targ.GetMass() << std::endl;
		std::cout << " Recoil mass (MeV/c^2)   = " << recoil.GetMass() << std::endl;
		std::cout << " Ejectile mass (MeV/c^2) = " << ejectile.GetMass() << std::endl;
		std::cout << " g.s. Q-value (MeV)      = " << Qgs << std::endl;

		kind = new Kindeux();
		kind->Initialize(beam.GetA(), targ.GetA(), recoil.GetA(), ejectile.GetA(), Qgs, 1, recoilStates);

		std::cout << " Inverse Kinematics:";
		if(kind->GetInverseKinematics())
			std::cout << " Yes\n";
		else
			std::cout << " No\n";
	}

	pack.Open("mcarlo.root", WriteReaction);

	std::cout << std::endl;

	if(Nejectile > 0){		
		Nwanted = (unsigned int)proper_value("Enter number of ejectile MC events: ", 0.0, true);

		// Process ejectile detectors
		if(Nwanted > 0){
			std::cout << "  Performing Monte Carlo test on ejectile detectors...\n";
			total_found = TestDetSetup(&pack, detectors, Nwanted, WriteReaction, beamspot, targangle, true, kind);
		
			std::cout << "  Found " << Nwanted << " ejectile events in " << total_found << " trials (" << 100.0*Nwanted/total_found << "%)\n\n";

			std::stringstream stream; stream << Nwanted;
			TNamed n1("EjectDet", stream.str().c_str());
			stream.str(""); stream << total_found;
			TNamed n2("EjectTot", stream.str().c_str());
			stream.str(""); stream << 100.0*Nwanted/total_found << " %";
			TNamed n3("EjectEff", stream.str().c_str());

			pack.file->cd();
			n1.Write();
			n2.Write();
			n3.Write();
		}
	}
	else{ std::cout << " Note: Found no ejectile detectors in detector file.\n\n"; }

	if(Nrecoil > 0){
		Nwanted = (unsigned int)proper_value("Enter number of recoil MC events: ", 0.0, true);

		// Process recoil detectors
		if(Nwanted > 0){
			std::cout << "  Performing Monte Carlo test on recoil detectors...\n";
			total_found = TestDetSetup(&pack, detectors, Nwanted, WriteReaction, beamspot, targangle, false, kind);
			
			std::cout << "  Found " << Nwanted << " recoil events in " << total_found << " trials (" << 100.0*Nwanted/total_found << "%)\n\n";

			std::stringstream stream; stream << Nwanted;
			TNamed n1("RecoilDet", stream.str().c_str());
			stream.str(""); stream << total_found;
			TNamed n2("RecoilTot", stream.str().c_str());
			stream.str(""); stream << 100.0*Nwanted/total_found << " %";
			TNamed n3("RecoilEff", stream.str().c_str());

			pack.file->cd();
			n1.Write();
			n2.Write();
			n3.Write();
		}
	}
	else{ std::cout << " Note: Found no recoil detectors in detector file.\n\n"; }

	std::cout << " Finished geometric efficiency test on detector setup...\n";
	
	for(std::vector<Primitive*>::iterator iter = detectors.begin(); iter != detectors.end(); iter++){
		delete *iter;
	}
	detectors.clear();

	return 0;
}
