// vikar.h
// Converted by FortranConvert v0.1
// Wed Feb 12 19:55:17 2014

#include <fstream>
#include <iostream>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"

#include "vikar_core.h"
#include "kindeux.h"
#include "materials.h"
#include "detectors.h"
#include "Structures.h"

#define VERSION "1.25c"

struct debugData{
	double var1, var2, var3;
	
	void Set(double v1, double v2, double v3){
		var1 = v1; var2 = v2; var3 = v3;
	}
};

template <typename T>
void SetName(std::vector<TNamed*> &named, std::string name_, const T &value_, std::string units_=""){
	std::stringstream stream;
	stream << value_;
	stream << " " << units_;
	
	named.push_back(new TNamed(name_.c_str(), stream.str().c_str()));
}

double MeV2MeVee(const double &Tp_){
	// Coefficients for NE102 (BC408)
	const double Acoeff = 0.95;
	const double Bcoeff = 8.00;
	const double Ccoeff = 0.10;
	const double Dcoeff = 0.90;
	
	return (Acoeff*Tp_ - Bcoeff*(1 - std::exp(-Ccoeff*std::pow(Tp_, Dcoeff))));
}

int main(int argc, char* argv[]){ 
	// A Monte-Carlo charged particle experiment simulation program - details below
	//
	// vikar 1.0 written by S.D. Pain on some date in 2004
	//
	// vikar 2.0 updated by S.D. Pain on 5/5/2006
	//   - Updated to employ non-isotropic angular distributions
	//
	// vikar 3.0 last updated by S.D. Pain on 6/13/2010
	//   - Major structural reworking
	//   - Updated to employ charged particle processing subroutine
	//   - Fixed to set detector properties (RadLength, conv) for SRIM tables
	//   - Fixed ranges array size in SRIM tables and main to match
	//
	// vikar 3.1 last updated by S.D. Pain on 10/15/2013
	//   - Added beam energy spread (SRIM input unchecked for functionality)
	//   - Added beamspot size parameter in input file
	//   - Applied beamspot size variation to annular and cylindrical detectors (1st-order application)
	//   - Applied version check to input file read/write
	//
	// vikar 3.11 last updated by S.D. Pain on 12/11/2013
	//   - resolution_cyl updated (v1.0 to v1.1) to improve beam-spot-size effects
	//   - product_proc updated (v1.0 to v1.1) to improve beam-spot-size effects for annular detectors
	//
	// Outstanding things for future versions
	//----------------------------------------
	// Write relativistic conversion routines
	// Put in recoil breakup
	// Add energy straggling (see NIMS 117, 125 (1974))
	// Tilt annular detectors
	// Add in gamma-rays
	// Make strips selectively resistive/non-resistive
	// Allow different strips to be hit on dE and E detectors
	// Allow ejectile excitations
	// Put in user defined energy for position resolution point
	// Add user selected detector material
	// Add different detector materials
	// Fix stuck problems for detecting both ejectiles and recoils when no coincidence
	// is required
	
	// Seed randomizer
	srand(time(NULL));
	
	// Main objects
	Kindeux kind; // Main kinematics object
	Target targ; // The physical target
	Efficiency bar_eff; // VANDLE bar efficiencies
	std::vector<Primitive*> vandle_bars; // Vector of Primitive detectors

	RangeTable beam_targ; // Range table for beam in target
	RangeTable eject_targ; // Pointer to the range table for ejectile in target
	RangeTable recoil_targ; // Pointer to the range table for recoil in target
	RangeTable *eject_tables = NULL; // Array of range tables for ejectile in various materials
	RangeTable *recoil_tables = NULL; // Array of range tables for recoil in various materials

	Particle recoil_part; // Recoil particle
	Particle eject_part;// Ejectile particle
	Particle beam_part; // Beam particle

	// Hit coordinates on the surface of a detector
	double hit_x, hit_y, hit_z;

	// Vectors and rotation matrices
	Vector3 ZeroVector; // The zero vector
	Vector3 Ejectile, Recoil, Gamma;
	Vector3 HitDetect1, HitDetect2;
	Vector3 EjectSphere, RecoilSphere;
	Vector3 lab_beam_focus; // The focus point for the beam. Non-cylindrical beam particles will originate from this point
	Vector3 lab_beam_start; // The originating point of the beam particle in the lab frame
	Vector3 lab_beam_trajectory; // The original trajectory of the beam particle before it enters the target
	Vector3 lab_beam_interaction; // The position of the reaction inside the target
	Vector3 lab_beam_stragtraject; // The angular straggled trajectory of the beam particle just before the reaction occurs
	Vector3 targ_surface; // The intersection point between the beam particle and the target surface (wrt beam focus)
	Vector3 interaction; // The interaction point inside the target (wrt beam focus)
	Matrix3 rotation_matrix; // The rotation matrix used to transform vectors from the beam particle frame to the lab frame

	double tof; // Neutron time of flight (s)
	double Zdepth; // Interaction depth inside of the target (m)
	double range_beam;
		
	unsigned int num_materials = 0;
	Material *materials = NULL; // Array of materials
	
	unsigned int targ_mat_id; // The ID number of the target material
	std::string targ_mat_name; // The name of the target material
	
	// Physics Variables
	unsigned int NRecoilStates = 0;
	std::vector<std::string> AngDist_fname; 
	double *ExRecoilStates = NULL;
	double *totXsect = NULL; 
	double gsQvalue = 0.0;

	// Energy variables
	double Ebeam = 0.0, Ebeam0 = 0.0;
	double ErecoilMod = 0.0, EejectMod = 0.0;
		
	// Beam variables
	double beamspot = 0.0; // Beamspot diameter (m) (on the surface of the target)
	double beamEspread = 0.0; // Beam energy spread (MeV)
	double beamAngdiv = 0.0; // Beam angular divergence (radians)

	double timeRes = 2E-9; // Pixie-16 time resolution (s)
	double BeamRate = 0.0; // Beam rate (1/s)

	// Background variables
	unsigned int backgroundRate = 0;
	unsigned int backgroundWait = 0;
	double detWindow = 0.0;
	bool bgPerDetection = false;
	
	std::string det_fname; // Detector filename
	std::string output_fname = "VIKAR.root";
	
	// Input/output variables
	unsigned int Ndetected = 0; // Total number of particles detected in VANDLE
	unsigned int Nwanted = 0; // Number of desired detections
	unsigned int Nsimulated = 0; // Total number of simulated particles
	unsigned int NdetHit = 0; // Total number of particles which collided with a bar
	unsigned int Nreactions = 0; // Total number of particles which react with the target
	int Ndet = 0; // Total number of detectors
	int NdetEject = 0; // Total number of ejectile detectors
	int NdetRecoil = 0; // Total number of recoil detectors
	clock_t timer; // Clock object for calculating time taken and remaining

	// Default options
	bool InverseKinematics = true;
	bool InCoincidence = true;
	bool WriteReaction = false;
	bool WriteDebug = false;
	bool PerfectDet = true;
	bool SupplyRates = false;
	bool BeamFocus = false;
	bool DoRutherford = false;
	bool CylindricalBeam = true;
	int ADists = 0;

	//------------------------------------------------------------------------
	//
	// Start of user input code
	//
	//------------------------------------------------------------------------

	std::cout << "\n####       #### ######## ####     ###         ###       #########      ########\n"; 
	std::cout << " ##         ##     ##     ##     ##          ## ##       ##     ##    ##      ##\n";
	std::cout << "  ##       ##      ##     ##    ##          ##   ##      ##     ##             ##\n";
	std::cout << "  ##       ##      ##     ##   ##           ##   ##      ##   ##              ##\n";
	std::cout << "   ##     ##       ##     #####            #########     #####               ##\n";
	std::cout << "   ##     ##       ##     ##   ##          ##     ##     ##   ##           ##\n";
	std::cout << "    ##   ##        ##     ##    ##        ##       ##    ##    ##         ##\n";
	std::cout << "    ##   ##        ##     ##     ##       ##       ##    ##     ##      ##\n";
	std::cout << "     ## ##         ##     ##      ##     ##         ##   ##      ##    ##\n";
	std::cout << "      ###       ######## ####      ###  ####       #### ####      ### ###########\n";

	std::cout << "\n VIKAR v " << VERSION << "\n"; 
	std::cout << " ==  ==  ==  ==  == \n\n"; 
	
	std::cout << " Welcome to NewVIKAR, the Virtual Instrumentation for Kinematics\n"; 
	std::cout << " And Reactions program, optimized for use with VANDLE bars\n\n";

	std::cout << " How about a nice cup of tea?\n"; 
	std::cout << " No? Well let me just load the input file then\n";
	std::cout << "\n ==  ==  ==  ==  == \n";
	sleep(1);

	std::ifstream input_file;
	if(argc >= 2){
		// Read an input file
		input_file.open(argv[1]);
		if(!input_file.good()){
			std::cout << " Error: Problem loading the input file\n";
			return 1;
		}
		
		// Specify the name of the output files
		if(argc >= 3){ output_fname = std::string(argv[2]); }
		
		unsigned int count = 0;
		std::string input;
		std::cout << "\n Reading from file " << argv[1] << std::endl;
		while(true){
			getline(input_file, input); input = Parse(input);
			if(input_file.eof()){ break; }
			if(input == ""){ continue; }
			
			if(count == 0){ 
				std::cout << "  Version: " << input << std::endl;
			}
			else if(count == 1){ 
				beam_part.SetZ(atof(input.c_str()));  
				std::cout << "  Beam-Z: " << beam_part.GetZ() << std::endl;
			}
			else if(count == 2){ 
				beam_part.SetA(atof(input.c_str()));
				std::cout << "  Beam-A: " << beam_part.GetA() << std::endl;
			}
			else if(count == 3){
				beam_part.SetMass(atof(input.c_str()));
				std::cout << "  Beam Mass: " << beam_part.GetMassAMU() << " amu\n";
			}
			else if(count == 4){ 
				targ.SetZ((double)atof(input.c_str()));
				std::cout << "  Target-Z: " << targ.GetZ() << std::endl;
			}
			else if(count == 5){ 
				targ.SetA((double)atof(input.c_str()));
				std::cout << "  Target-A: " << targ.GetA() << std::endl;
			}
			else if(count == 6){
				targ.SetMass(atof(input.c_str()));
				std::cout << "  Target Mass: " << targ.GetMassAMU() << " amu\n";
			}
			else if(count == 7){ 
				recoil_part.SetZ((double)atof(input.c_str()));
				std::cout << "  Recoil-Z: " << recoil_part.GetZ() << std::endl;
			}
			else if(count == 8){ 
				recoil_part.SetA((double)atof(input.c_str()));
				std::cout << "  Recoil-A: " << recoil_part.GetA() << std::endl;
			}
			else if(count == 9){
				recoil_part.SetMass(atof(input.c_str()));
				std::cout << "  Recoil Mass: " << recoil_part.GetMassAMU() << " amu\n";
			}
			else if(count == 10){ 
				eject_part.SetZ(atof(input.c_str()));
				recoil_part.SetZ(beam_part.GetZ() + targ.GetZ() - eject_part.GetZ());
				std::cout << "  Ejectile-Z: " << eject_part.GetZ() << std::endl;
			}
			else if(count == 11){ 
				eject_part.SetA(atof(input.c_str())); 
				std::cout << "  Ejectile-A: " << eject_part.GetA() << std::endl;
			}
			else if(count == 12){
				eject_part.SetMass(atof(input.c_str()));
				std::cout << "  Ejectile Mass: " << eject_part.GetMassAMU() << " amu\n";
				
				// Set the reaction Q-value.
				gsQvalue = (beam_part.GetMass()+targ.GetMass())-(recoil_part.GetMass()+eject_part.GetMass()); 
				std::cout << "  G.S. Q-Value: " << gsQvalue << " MeV\n";
				
				// Set inverse or normal kinematics.
				if(beam_part.GetA() > targ.GetA()){ 
					InverseKinematics = true; 
					std::cout << "  Inverse Kinematics: Yes\n";
				}
				else{ 
					InverseKinematics = false; 
					std::cout << "  Inverse Kinematics: No\n";
				}
			}
			else if(count == 13){ 
				Ebeam0 = atof(input.c_str());  
				std::cout << "  Beam Energy: " << Ebeam0 << " MeV\n";
			}
			else if(count == 14){ 
				SetBool(input, "  Cylindrical Beam", CylindricalBeam);
				getline(input_file, input); input = Parse(input);
				beamspot = atof(input.c_str());
				if(CylindricalBeam){ std::cout << "  Beam Spot Diameter: " << beamspot << " mm\n"; }
				else{ std::cout << "  Beam Spot FWHM: " << beamspot << " mm\n"; }
				beamspot = beamspot/1000.0; // in meters
			}
			else if(count == 15){ 
				beamAngdiv = atof(input.c_str());  
				std::cout << "  Beam Angular Divergence: " << beamAngdiv << " degrees\n";
				beamAngdiv *= deg2rad; // in radians
			}
			else if(count == 16){ 
				beamEspread = atof(input.c_str());  
				std::cout << "  Beam Energy Spread: " << beamEspread << " MeV\n";
			}
			else if(count == 17){ 
				// Recoil excited state information
				NRecoilStates = atoi(input.c_str()) + 1;
				std::cout << "  No. Excited States: " << NRecoilStates-1 << std::endl;
				ExRecoilStates = new double[NRecoilStates]; ExRecoilStates[0] = 0.0;
				totXsect = new double[NRecoilStates]; totXsect[0] = 0.0;
				std::cout << "   Recoil Ground State: 0.0 MeV\n";
				for(unsigned int i = 1; i < NRecoilStates; i++){
					getline(input_file, input); input = Parse(input);
					ExRecoilStates[i] = atof(input.c_str());
					std::cout << "   Recoil Excited State " << i << ": " << ExRecoilStates[i] << " MeV\n";
					totXsect[i] = 0.0;
				}
			}
			else if(count == 18){ 
				// Angular distribution information
				ADists = atoi(input.c_str());
				if(ADists == 1){ // Read in filename for each state's angular distribution.
					std::cout << "  Supply Angular Distributions: Yes\n";
					
					// Read the filenames.
					for(unsigned int i = 0; i < NRecoilStates; i++){
						getline(input_file, input); input = Parse(input);
						if(!DoRutherford && input != "RUTHERFORD"){
							AngDist_fname.push_back("./xsections/" + input);
							if(i == 0){ std::cout << "   Distribution for ground state: " << AngDist_fname[i] << std::endl; }
							else{ std::cout << "   Distribution for state " << i+1 << ": " << AngDist_fname[i] << std::endl; }
						}
						else{
							DoRutherford = true;
							std::cout << "   Using Rutherford scattering\n";
						}
					}
				
					// Supply beam rate information.
					getline(input_file, input); input = Parse(input);
					if(SetBool(input, "  Calculate Rates", SupplyRates)){ 
						getline(input_file, input); input = Parse(input);
						BeamRate = atof(input.c_str());
						std::cout << "   Beam Rate: " << BeamRate << " pps\n";
					}
				}
				else if(ADists == 2){
					std::cout << "  Supply Angular Distributions: Yes\n";
					
					// Supply the rates relative to the ground state. i.e. a value of
					// 1.0 will generate 1 particle for each ground state particle. A
					// value of 2.0 will generate 2 particles for each ground state and so on.
					for(unsigned int i = 0; i < NRecoilStates; i++){
						getline(input_file, input); input = Parse(input);
						AngDist_fname.push_back(input);
						if(i == 0){ std::cout << "   Production rate ground state: " << input << " per event.\n"; }
						else{ std::cout << "   Production rate for state " << i+1 << ": " << input << " per event.\n"; }
					}
				}
				else{ std::cout << "  Supply Angular Distributions: No\n"; }
			}
			else if(count == 19){
				// Target material
				targ_mat_name = input;
			}
			else if(count == 20){ 
				// Target thickness
				if(targ_mat_name != "NONE"){ // Normal material.
					targ.SetThickness((double)atof(input.c_str()));
					std::cout << "  Target Thickness: " << targ.GetThickness() << " mg/cm^2\n";	
				}
				else{ // Target material energy loss disabled.
					targ.SetRealThickness((double)atof(input.c_str()));
					std::cout << "  Target Thickness: " << targ.GetRealThickness() << " m\n";	
				}		
			}
			else if(count == 21){ 
				// Target angle wrt beam axis
				targ.SetAngle((double)atof(input.c_str())*deg2rad);
				std::cout << "  Target Angle: " << targ.GetAngle()*rad2deg << " degrees\n";
			}
			else if(count == 22){ 
				// Load the small, medium, and large bar efficiencies
				// Efficiency index 0 is the underflow efficiency (for energies below E[0])
				// Efficiency index N is the overflow efficiency (for energies greater than E[N])
				if(!SetBool(input, "  Perfect Detector", PerfectDet)){ 
					// Load small bar efficiency data
					getline(input_file, input); input = "./efficiency/" + Parse(input); 
					std::cout << "   Found " << bar_eff.ReadSmall(input.c_str()) << " small bar data points in file " << input << "\n";
					
					// Load medium bar efficiency data
					getline(input_file, input); input = "./efficiency/" + Parse(input); 
					std::cout << "   Found " << bar_eff.ReadMedium(input.c_str()) << " medium bar data points in file " << input << "\n";
					
					// Load large bar efficiency data
					getline(input_file, input); input = "./efficiency/" + Parse(input); 
					std::cout << "   Found " << bar_eff.ReadLarge(input.c_str()) << " large bar data points in file " << input << "\n";
				}
			}
			else if(count == 23){ 
				// Load detector setup from a file
				det_fname = "./detectors/" + input;
				std::cout << "  Detector Setup Filename: " << det_fname << std::endl;
			}
			else if(count == 24){ 
				// Desired number of detections
				Nwanted = atol(input.c_str());
				std::cout << "  Desired Detections: " << Nwanted << std::endl; 
			}
			else if(count == 25){
				// Background rate (per detection event)
				backgroundRate = atol(input.c_str());
				if(backgroundRate > 0){ // Get the detection ToF window
					getline(input_file, input); input = Parse(input);
					detWindow = atof(input.c_str());
					backgroundWait = backgroundRate;
					getline(input_file, input); input = Parse(input);
					SetBool(input, bgPerDetection);
					if(bgPerDetection){	std::cout << "  Background Rate: " << backgroundRate << " events per detection\n"; }
					else{ std::cout << "  Background Rate: " << backgroundRate << " events per recoil\n"; }
					std::cout << "  Background Time Width: " << detWindow << " ns\n";
				}
				else{ std::cout << "  Background Rate: NONE\n"; }
			}
			else if(count == 26){
				// Require ejectile and recoil particle coincidence?
				SetBool(input, "  Require Particle Coincidence", InCoincidence);
			}
			else if(count == 27){
				// Write Reaction data to file?
				SetBool(input, "  Write Reaction Info", WriteReaction);
			}
			else if(count == 28){
				// Write Debug data to file?
				SetBool(input, "  Write Debug Info", WriteDebug);
			}
			
			count++;
		}
		
		input_file.close();
		if(count <= 28){ std::cout << " Warning! The input file is invalid. Check to make sure input is correct\n"; }
	}
	else{
		std::cout << " FATAL ERROR! Missing required variable! Aborting...\n";
		return 1;
	}
		
	if(beam_part.GetA()+targ.GetA() != recoil_part.GetA()+eject_part.GetA()){
		std::cout << "\n FATAL ERROR! Mass number is NOT conserved! Aborting...\n";
		return 1;
	}
		
	std::cout << "\n ==  ==  ==  ==  == \n\n";

	// Make sure the input variables are correct
	if(!Prompt(" Are the above settings correct?")){
		std::cout << "  ABORTING...\n";
		return 1;
	}

	std::cout << "\n Initializing main simulation Kindeux object...\n";

	// Initialize kinematics object
	kind.Initialize(beam_part.GetA(), targ.GetA(), recoil_part.GetA(), eject_part.GetA(), gsQvalue, NRecoilStates, ExRecoilStates);

	// Read the detector setup file
	std::cout << " Reading in NewVIKAR detector setup file...\n";
	Ndet = ReadDetFile(det_fname.c_str(), vandle_bars);
	if(Ndet < 0){ // Failed to load setup file
		std::cout << " Error: failed to load detector setup file!\n";
		return 1; 
	}
	else if(Ndet == 0){ std::cout << " Error: Found no detectors in the detector setup file!\n"; } // Check there's at least 1 detector!
	
	std::vector<std::string> needed_materials;
	for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
		if(!IsInVector((*iter)->GetMaterialName(), needed_materials)){
			needed_materials.push_back((*iter)->GetMaterialName());
		}
		
		if((*iter)->IsEjectileDet()){ NdetEject++; }
		if((*iter)->IsRecoilDet()){ NdetRecoil++; }
	}

	// Report on how many detectors were read in
	std::cout << " Found " << NdetEject << " ejectile and " << NdetRecoil << " recoil detectors in file " << det_fname << std::endl;

	bool use_target_eloss = true;

	// Load VIKAR material files
	targ_mat_id = 0;
	std::ifstream material_names("./materials/names.in");
	if(material_names.good()){
		std::cout << "\n Loading VIKAR material files...\n";
		std::vector<std::string> names;
		std::string line;
		while(true){
			getline(material_names, line);
			if(material_names.eof()){ break; }
			if(line[0] == '#'){ continue; } // Commented line
			
			line = Parse(line);
			names.push_back(line);
		}
		
		materials = new Material[names.size()+1];
		if(eject_part.GetZ() > 0){ eject_tables = new RangeTable[names.size()+1]; }
		if(recoil_part.GetZ() > 0){ recoil_tables = new RangeTable[names.size()+1]; }
		num_materials = 1; // Default CD2 material
		for(std::vector<std::string>::iterator iter = names.begin(); iter != names.end(); iter++){
			materials[num_materials].ReadMatFile(iter->c_str());
			num_materials++;
		}
		
		std::cout << " Successfully loaded " << num_materials-1 << " materials\n";
		if(targ_mat_name == "CD2"){ } // Use the default target
		else if(targ_mat_name == "NONE"){ 
			std::cout << "  Target Material: DISABLED\n";
			use_target_eloss = false;
		}
		else{
			std::cout << "  Target Material: " << targ_mat_name;
			for(unsigned int i = 1; i < num_materials; i++){
				if(targ_mat_name == materials[i].GetName()){
					targ_mat_id = i;
					break;
				}
			}
			if(targ_mat_id == 0){ 
				std::cout << " (not found)"; 
				use_target_eloss = false;			
			}
			std::cout << std::endl;
		}
	}
	else{ 
		std::cout << " Warning! Failed to load the file ./materials/names.in\n"; 
		materials = new Material[1];
		if(eject_part.GetZ() > 0){ eject_tables = new RangeTable[1]; }
		if(recoil_part.GetZ() > 0){ recoil_tables = new RangeTable[1]; }
	}

	// Setup default CD2 material
	materials[0].SetName("CD2");
	materials[0].Init(2);
	materials[0].SetDensity(1.06300);
	materials[0].SetMolarMass(16.038904);

	unsigned int num_per_molecule[2] = {1, 2};
	double element_Z[2] = {6, 1};
	double element_A[2] = {12, 2};
	materials[0].SetElements(num_per_molecule, element_Z, element_A);
	
	if(use_target_eloss){
		if(targ_mat_id == 0){
			std::cout << "  Target Material: CD2\n";
		}
	
		targ.SetDensity(materials[targ_mat_id].GetDensity());
		targ.SetRadLength(materials[targ_mat_id].GetRadLength());
		std::cout << "  Target Radiation Length: " << targ.GetRadLength() << " mg/cm^2\n\n";
	
		// Calculate the stopping power table for the reactino particles in the target
		if(beam_part.GetZ() > 0){ // The beam is a charged particle (not a neutron)
			std::cout << " Calculating range table for beam in " << materials[targ_mat_id].GetName() << "...";
			beam_targ.Init(100, 0.1, (Ebeam0+2*beamEspread), beam_part.GetZ(), beam_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
		if(eject_part.GetZ() > 0){
			std::cout << " Calculating range table for ejectile in " << materials[targ_mat_id].GetName() << "...";
			eject_targ.Init(100, 0.1, (Ebeam0+2*beamEspread), eject_part.GetZ(), eject_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
		if(recoil_part.GetZ() > 0){
			std::cout << " Calculating range table for recoil in " << materials[targ_mat_id].GetName() << "...";
			recoil_targ.Init(100, 0.1, (Ebeam0+2*beamEspread), recoil_part.GetZ(), recoil_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
	}

	// Calculate the stopping power table for the ejectiles in the materials
	if(eject_part.GetZ() > 0){ // The ejectile is a charged particle (not a neutron)
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating ejectile range table for " << materials[i].GetName() << "...";
			eject_tables[i].Init(100, eject_part.GetKEfromV(0.02*c), (Ebeam0+2*beamEspread), eject_part.GetZ(), eject_part.GetA()/mev2amu, &materials[i]);
			std::cout << " Done!\n";
		}
	}


	// Calculate the stopping power table for the recoils in the materials
	if(recoil_part.GetZ() > 0){ // The recoil is a charged particle (not a neutron)
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating recoil range table for " << materials[i].GetName() << "...";
			recoil_tables[i].Init(100, recoil_part.GetKEfromV(0.02*c), (Ebeam0+2*beamEspread), recoil_part.GetZ(), recoil_part.GetA()/mev2amu, &materials[i]);
			std::cout << " Done!\n";
		}
	}
	
	// Calculate the beam focal point (if it exists)
	if(beamAngdiv >= 0.000174532925199){ // Beam focus is upstream of target
		lab_beam_focus = Vector3(0.0, 0.0, -(beamspot/(2.0*std::tan(beamAngdiv))+targ.GetRealZthickness()/2.0));
		std::cout << " Beam focal point at Z = " << lab_beam_focus.axis[2] << " m\n";
		BeamFocus = true;
	}
	if(beamAngdiv <= -0.000174532925199){ // Beam focus is downstream of target
		lab_beam_focus = Vector3(0.0, 0.0, (beamspot/(2.0*std::tan(beamAngdiv))+targ.GetRealZthickness()/2.0));
		std::cout << " Beam focal point at Z = " << lab_beam_focus.axis[2] << " m\n";
		BeamFocus = true;
	}

	// For cylindrical beams, the beam direction is given by the z-axis
	if(!BeamFocus){ lab_beam_trajectory = Vector3(0.0, 0.0, 1.0); }

	std::cout << "\n Setting detector material types...\n";
	for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){ // Set the detector material for energy loss calculations
		for(unsigned int j = 0; j < num_materials; j++){
			if((*iter)->GetMaterialName() == materials[j].GetName()){
				if(((*iter)->IsRecoilDet() && recoil_part.GetZ() > 0) || ((*iter)->IsEjectileDet() && eject_part.GetZ() > 0)){ 
					// Only set detector to use material if the particle it is responsible for detecting has
					// a Z greater than zero. Particles with Z == 0 will not have calculated range tables
					// and thus cannot use energy loss considerations.
					(*iter)->SetMaterial(j);  
				}
				break; 
			}
		}
	}

	//std::cout << "\n ==  ==  ==  ==  == \n\n";

	if((ADists == 1 || ADists == 2) && !DoRutherford){ 
		if(ADists == 1){
			std::cout << "\n Loading state angular distribution files...\n";
			if(kind.SetDist(AngDist_fname, materials[targ_mat_id].GetTotalElements(), materials[targ_mat_id].GetDensity(), BeamRate)){
				// Successfully set the angular distributions
				std::cout << " Successfully loaded angular distributions.\n";
				kind.Print();
			}
			else{
				std::cout << "  Warning! Could not properly initialize distributions.\n";
				std::cout << "  Note: Setting all energy states to isotropic!\n";
			}
		}
		else{
			std::cout << "\n Setting relative intensities for recoil states...\n";
			if(kind.SetDist(AngDist_fname)){
				// Successfully set the angular distributions
				std::cout << " Successfully set relative intensities for recoil states.\n";
			}
			else{
				std::cout << "  Warning! Could not properly initialize distributions.\n";
				std::cout << "  Note: Setting all energy states to isotropic!\n";
			}
		}
	}
	else if(DoRutherford){
		double e = 1.60217657E-19; // C
		double k = 8.987551E9; // N*m^2/C^2
		double coefficient = k*beam_part.GetZ()*targ.GetZ()*e*e/(4.0*Ebeam0*1.60218E-13); // m^2
		coefficient = coefficient * coefficient;
		std::cout << "\n Generating Rutherford distribution...\n";
		std::cout << "  Z-1 = " << beam_part.GetZ()*e << " C\n";
		std::cout << "  Z-2 = " << targ.GetZ()*e << " C\n";
		std::cout << "  Energy = " << Ebeam0*1.60218E-13 << " J\n";
		std::cout << "  Coefficient: " << coefficient << " m^2\n";
		
		if(NRecoilStates > 1){
			std::cout << "   Warning! Cannot set distribution to Rutherford for excited states.\n";
			std::cout << "   Note: Setting number of excited states to zero!\n";
			NRecoilStates = 1;
		}
		kind.SetRutherford(coefficient);
	}

	std::cout << "\n ==  ==  ==  ==  == \n\n";

	// Last chance to abort
	if(!Prompt(" Setup is complete. Is everything correct?")){
		std::cout << "  ABORTING...\n";
		return 1;
	}

	//---------------------------------------------------------------------------
	// End of Input Section
	//---------------------------------------------------------------------------
		
	// Root stuff
	TFile *file = new TFile(output_fname.c_str(), "RECREATE");
	TTree *VIKARtree = new TTree("VIKAR", "VIKAR output tree");
	TTree *DEBUGtree = NULL;
	
	EjectObjectStructure EJECTdata;
	RecoilObjectStructure RECOILdata;
	ReactionObjectStructure REACTIONdata;
	debugData DEBUGdata;
	
	VIKARtree->Branch("Eject", &EJECTdata);
	VIKARtree->Branch("Recoil", &RECOILdata);
	if(WriteReaction){
		VIKARtree->Branch("Reaction", &REACTIONdata);
	}
	if(WriteDebug){ 
		DEBUGtree = new TTree("DEBUG", "VIKAR debug tree");
		DEBUGtree->Branch("Debug", &DEBUGdata, "var1/D:var2/D:var3/D"); 
	}

	// Write reaction info to the file.
	std::vector<TNamed*> named;
	SetName(named, "Version", VERSION);
	SetName(named, "Beam-Z", beam_part.GetZ());
	SetName(named, "Beam-A", beam_part.GetA());
	SetName(named, "Beam Mass", beam_part.GetMassAMU(), "amu");
	SetName(named, "Target-Z", targ.GetZ());
	SetName(named, "Target-A", targ.GetA());
	SetName(named, "Target Mass", targ.GetMassAMU(), "amu");
	SetName(named, "Recoil-Z", recoil_part.GetZ());
	SetName(named, "Recoil-A", recoil_part.GetA());
	SetName(named, "Recoil Mass", recoil_part.GetMassAMU(), "amu");
	SetName(named, "Ejectile-Z", eject_part.GetZ());
	SetName(named, "Ejectile-A", eject_part.GetA());
	SetName(named, "Ejectile Mass", eject_part.GetMassAMU(), "amu");
	if(InverseKinematics){ SetName(named, "Inverse Kinematics?", "Yes"); }
	else{ SetName(named, "Inverse Kinematics?", "No"); }
	SetName(named, "Beam Energy", Ebeam0, "MeV");	
	if(CylindricalBeam){
		SetName(named, "Gaussian Beam", "No");	
		SetName(named, "Beamspot Dia.", beamspot, "m");	
	}
	else{
		SetName(named, "Gaussian Beam", "Yes");	
		SetName(named, "Beamspot FWHM", beamspot, "m");	
	}
	SetName(named, "Beam Divergence", beamAngdiv*rad2deg, "deg");
	SetName(named, "Beam dE", beamEspread, "MeV");
	SetName(named, "G.S. Q-Value", gsQvalue, "MeV");
	SetName(named, "No. Excited", NRecoilStates-1, "MeV");
	SetName(named, "Recoil G.S.", "0.0", "MeV");
	for(unsigned int i = 1; i < NRecoilStates; i++){
		std::stringstream stream; stream << i;
		SetName(named, "Recoil Ex. State "+stream.str(), ExRecoilStates[i], "MeV");
	}
	if(ADists == 1){ 
		SetName(named, "Supply Distributions?", "Yes, 1"); 
		for(unsigned int i = 0; i < NRecoilStates; i++){
			std::stringstream stream; stream << i;
			SetName(named, "State "+stream.str()+" Dist.", AngDist_fname[i]);
		}
	}
	else if(ADists == 2){
		SetName(named, "Supply Distributions?", "Yes, 2");
		for(unsigned int i = 0; i < NRecoilStates; i++){
			std::stringstream stream; stream << i;
			SetName(named, "State "+stream.str()+" Rate", AngDist_fname[i], "per event");
		}
	}
	else{ SetName(named, "Supply Distributions?", "No"); }
	SetName(named, "Target Material", targ_mat_name);
	if(targ_mat_name != "NONE"){ SetName(named, "Target Thickness", targ.GetThickness(), "mg/cm^2"); }	
	else{ SetName(named, "Target Thickness", targ.GetRealThickness(), "m"); }	
	SetName(named, "Target Angle", targ.GetAngle()*rad2deg, "deg");
	if(PerfectDet){ SetName(named, "Perfect Detectors?", "Yes"); }
	else{ SetName(named, "Perfect Detectors?", "No"); }
	SetName(named, "Detector Filename", det_fname);
	SetName(named, "Number Detections", Nwanted);
	if(backgroundRate > 0){ 
		if(bgPerDetection){ SetName(named, "Background Rate", backgroundRate, "per detection"); }
		else{ SetName(named, "Background Rate", backgroundRate, "per recoil"); }
		SetName(named, "Background Window", detWindow, "ns");
	}
	else{ SetName(named, "Background Rate", "NONE"); }
	if(InCoincidence){ SetName(named, "Recoil Coincidence?", "Yes"); }
	else{ SetName(named, "Recoil Coincidence?", "No"); }
	if(WriteReaction){ SetName(named, "Write Reaction?", "Yes"); }
	else{ SetName(named, "Write Reaction?", "No"); }
	if(WriteDebug){ SetName(named, "Write Debug?", "Yes"); }
	else{ SetName(named, "Write Debug?", "No"); }

	for(std::vector<TNamed*>::iterator iter = named.begin(); iter != named.end(); iter++){
		(*iter)->Write();
		delete (*iter);
	}
	named.clear();

	// Begin the simulation
	std::cout << " ---------- Simulation Setup Complete -----------\n"; 
	std::cout << "\n Beginning simulating " << Nwanted << " events....\n"; 

	//---------------------------------------------------------------------------
	// The Event Loop
	// ==  ==  ==  ==  ==  ==  == 
	// (Just to make it obvious)
	//---------------------------------------------------------------------------

	Vector3 temp_vector;
	Vector3 temp_vector_sphere;
	Vector3 dummy_vector;
	double dummy_t1, dummy_t2;
	double dist_traveled = 0.0, QDC = 0.0;
	double penetration = 0.0, fpath1 = 0.0, fpath2 = 0.0;
	float totTime = 0.0;
	char counter = 1;
	bool flag = false;
	bool hit = false;
	unsigned int chunk = Nwanted/10;
	unsigned int beam_stopped = 0;
	unsigned int recoil_stopped = 0;
	unsigned int eject_stopped = 0;
	timer = clock();
	
	bool proc_eject = false;
	
	// Struct for storing reaction information.
	reactData rdata;
	
	while(Ndetected < Nwanted){
		// ****************Time Estimate**************
		if(flag && (Ndetected % chunk == 0)){
			flag = false;
			totTime = (float)(clock()-timer)/CLOCKS_PER_SEC;
			std::cout << "\n ------------------------------------------------\n"; 
			std::cout << " Number of particles Simulated: " << Nsimulated << std::endl; 
			std::cout << " Number of particles Detected: " << Ndetected << std::endl; 
			if(SupplyRates && ADists){ std::cout << " Number of Reactions: " << Nreactions << std::endl; }
		
			std::cout << " " << Ndetected*100.0/Nwanted << "% of simulation complete...\n"; 
			if(PerfectDet){ std::cout << "  Detection Efficiency: " << Ndetected*100.0/Nreactions << "%\n"; }
			else{
				std::cout << "  Geometric Efficiency: " << NdetHit*100.0/Nreactions << "%\n";
				std::cout << "  Detection Efficiency: " << Ndetected*100.0/Nreactions << "%\n"; 
			}
			if(SupplyRates){ std::cout << "  Beam Time: " << Nsimulated/BeamRate << " seconds\n"; }
		
			std::cout << "  Simulation Time: " << totTime << " seconds\n";
			std::cout << "  Time reamining: " << (totTime/counter)*(10-counter) << " seconds\n";
			counter++; 
		}

		if(backgroundWait != 0){ // Simulating background events
			backgroundWait--;

			// Process the background event for each detector
			for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
				if(!((*iter)->IsEjectileDet() || (*iter)->IsRecoilDet())){ continue; } // This detector cannot detect particles

				// Select the "tof" of the background event
				tof = (double)frand(0, detWindow)*(1E-9);
			
				// Select a random point insde the detector
				(*iter)->GetRandomPointInside(temp_vector);
				Cart2Sphere(temp_vector_sphere);
			
				// Calculate the apparent energy of the particle using the tof
				if((*iter)->IsEjectileDet()){
					double dummyE = 0.5*kind.GetMejectMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					EJECTdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], temp_vector_sphere.axis[1]*rad2deg,
									 temp_vector_sphere.axis[2]*rad2deg, dummyE, tof*(1E9), 0.0, 0.0, 0.0, 0.0, (*iter)->GetLoc(), true);
					VIKARtree->Fill(); 
					EJECTdata.Zero();
				}
				else if((*iter)->IsRecoilDet()){
					double dummyE = 0.5*kind.GetMrecoilMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					RECOILdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], RecoilSphere.axis[1]*rad2deg,
									  RecoilSphere.axis[2]*rad2deg, dummyE, tof*(1E9), 0.0, 0.0, 0.0, 0.0, (*iter)->GetLoc(), true);
					VIKARtree->Fill();
					RECOILdata.Zero();
				}
			}
			
			continue;
		}
		else{ // A reaction occured
			if(!bgPerDetection){ backgroundWait = backgroundRate; }

			Nsimulated++; 
		
			// Simulate a beam particle before entering the target
			// Randomly select a point uniformly distributed on the beamspot
			// Calculate where the beam particle reacts inside the target
			// beamspot as well as the distance traversed through the target
			if(BeamFocus){ 
				// In this case, lab_beam_focus is the originating point of the beam particle
				// (or the terminating point for beams focused downstream of the target)
				// The direction is given by the cartesian vector 'lab_beam_trajectory'
				if(CylindricalBeam){ RandomCircle(beamspot, lab_beam_focus.axis[2], lab_beam_trajectory); } // Cylindrical beam
				else{ RandomGauss(beamspot, lab_beam_focus.axis[2], lab_beam_trajectory); } // Gaussian beam
				Zdepth = targ.GetInteractionDepth(lab_beam_focus, lab_beam_trajectory, targ_surface, lab_beam_interaction);
			}
			else{ 
				// In this case, lab_beam_start stores the originating point of the beam particle
				// The direction is given simply by the +z-axis
				// The 1m offset ensures the particle originates outside the target
				if(CylindricalBeam){ RandomCircle(beamspot, 1.0, lab_beam_start); } 
				else{ RandomGauss(beamspot, 1.0, lab_beam_start); } 
				lab_beam_start.axis[2] *= -1; // This is done to place the beam particle upstream of the target
				Zdepth = targ.GetInteractionDepth(lab_beam_start, lab_beam_trajectory, targ_surface, lab_beam_interaction);
			}	

			// Calculate the beam particle energy, varied with energy spread (in MeV)
			Ebeam = Ebeam0 + rndgauss0(beamEspread); 

			if(use_target_eloss){ // Calculate energy loss in the target
				// Calculate the beam particle range in the target (in m)
				range_beam = beam_targ.GetRange(Ebeam);
		
				// Calculate the new energy
				if(range_beam - Zdepth <= 0.0){ // The beam stops in the target (no reaction)
					if(beam_stopped == 10000){
						std::cout << "\n ATTENTION!\n";
						std::cout << "  A large number of beam particles (" << 100.0*beam_stopped/Nsimulated << "%) have stopped in the target!\n";
						std::cout << "  A high percentage of stopped particles could mean that the target is too thick.\n";
						std::cout << "  If this is the case, change the target thickness and restart the simulation.\n";
				
						std::cout << "\n ------------------------------------------------\n";
						std::cout << " Dumping target information!!!\n\n";
				
						materials[targ_mat_id].Print();
						std::cout << std::endl;
				
						std::cout << " Target thickness: " << targ.GetRealThickness() << " m\n";
						std::cout << " Effective thickness: " << targ.GetRealZthickness() << " m\n\n";
						
						std::cout << " Tracing ray through target... \n";
						std::cout << "  Target thickness: " << targ.GetPrimitive()->GetApparentThickness(Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, 1.0), dummy_vector, dummy_t1, dummy_t2) << " m\n";
						std::cout << "  Front face intersect = (" << dummy_vector.Dump() << ")\n";
						std::cout << "  Back face intersect = (" << Vector3(0.0, 0.0, -1 + dummy_t2).Dump() << ")\n";
					}
					beam_stopped++;
			
					continue; 
				}
				rdata.Ereact = beam_targ.GetEnergy(range_beam - Zdepth);
		
				// Determine the angle of the beam particle's trajectory at the
				// interaction point, due to angular straggling and the incident trajectory.
				targ.AngleStraggling(lab_beam_trajectory, beam_part.GetA(), beam_part.GetZ(), Ebeam, lab_beam_stragtraject);
			}
			else{ 
				rdata.Ereact = Ebeam;
				lab_beam_stragtraject = lab_beam_trajectory;
			}

			// the 2 body kinematics routine to generate the ejectile and recoil
			if(kind.FillVars(rdata, EjectSphere, RecoilSphere)){ Nreactions++; }
			else{ continue; } // A reaction did not occur

			// Convert the reaction vectors to cartesian coordinates
			// EjectSphere and RecoilSphere are unit vectors (no need to normalize)
			Sphere2Cart(EjectSphere, Ejectile); 
			Sphere2Cart(RecoilSphere, Recoil);
		
			// Transform the ejectile and recoil vectors (cartesian) from the beam trajectory frame into the Lab frame
			// This transformation will overwrite the Ejectile and Recoil vectors
			//rotation_matrix.SetRotationMatrixCart(lab_beam_trajectory); // Turn OFF angular straggling effects
			rotation_matrix.SetRotationMatrixCart(lab_beam_stragtraject); // Turn ON angular straggling effects
			rotation_matrix.Transform(Ejectile);
			rotation_matrix.Transform(Recoil);

			// Set the outgoing energy of the ejectile and recoil.
			EejectMod = rdata.Eeject;
			ErecoilMod = rdata.Erecoil;

			// Calculate the energy loss for the ejectile and recoil in the target.
			if(use_target_eloss){
				// Calculate the new energy of the ejectile.
				if(eject_part.GetZ() > 0){ 
					targ.GetPrimitive()->IntersectPrimitive(lab_beam_interaction, Ejectile, dummy_vector, Zdepth, dummy_t2);
					EejectMod = eject_targ.GetNewE(rdata.Eeject, Zdepth); 
				}
				
				// Calculate the new energy of the recoil.
				if(recoil_part.GetZ() > 0){ 
					targ.GetPrimitive()->IntersectPrimitive(lab_beam_interaction, Recoil, dummy_vector, Zdepth, dummy_t2);
					ErecoilMod = recoil_targ.GetNewE(rdata.Erecoil, Zdepth); 
				}
			}
		
			// Write debug information to the output file.
			if(WriteDebug){ 
				DEBUGdata.Set(Ejectile.axis[0], Ejectile.axis[1], Ejectile.axis[2]);
				DEBUGtree->Fill();
			}
		}

		// Process the reaction products
		for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
			if((*iter)->IsEjectileDet()){ // This detector can detect ejectiles
				if(EejectMod <= 0.0){ continue; } // There still may be recoil detectors left
				proc_eject = true; 
			}
			else if((*iter)->IsRecoilDet()){ // This detector can detect recoils
				if(ErecoilMod <= 0.0){ continue; } // There still may be ejectile detectors left
				proc_eject = false; 
			}
			else{ continue; } // This detector cannot detect particles, so skip it

process:			
			if(proc_eject){ 
				hit = (*iter)->IntersectPrimitive(lab_beam_interaction, Ejectile, HitDetect1, fpath1, fpath2); 
				
				// Calculate the vector pointing from the first intersection point to the second
				if(fpath1 >= 0.0 && fpath2 >= 0.0){  
					// The ray originates outside the detector. In this case, HitDetect1 represents the vector
					// from the origin to the point at which the ray intersects the detector surface.
					temp_vector = ((lab_beam_interaction + Ejectile*fpath2)-HitDetect1); 
				}
				else{ 
					// The ray originates within the detector. HitDetect1 already represents the vector from
					// the origin to the intersection of the surface of the detector.
					temp_vector = HitDetect1; 
				}
			}
			else{ 
				hit = (*iter)->IntersectPrimitive(lab_beam_interaction, Recoil, HitDetect1, fpath1, fpath2); 
				
				// Calculate the vector pointing from the first intersection point to the second
				if(fpath1 >= 0.0 && fpath2 >= 0.0){ 
					// The ray originates outside the detector. In this case, HitDetect1 represents the vector
					// from the origin to the point at which the ray intersects the detector surface.
					temp_vector = ((lab_beam_interaction + Recoil*fpath2)-HitDetect1); 
				}
				else{ 
					// The ray originates within the detector. HitDetect1 already represents the vector from
					// the origin to the intersection of the surface of the detector.
					temp_vector = HitDetect1; 
				}
			}
			
			// If a geometric hit was detected, process the particle
			if(hit){
				NdetHit++; 
				
				// The time of flight is the time it takes the particle to traverse the distance
				// from the target to the intersection point inside the detector
				penetration = temp_vector.Length(); // Total distance traveled through detector	

				// Solve for the energy deposited in the material.
				if((*iter)->UseMaterial()){ // Do energy loss and range considerations
					if(proc_eject){
						if(eject_part.GetZ() > 0){ // Calculate energy loss for the ejectile in the detector
							QDC = EejectMod - eject_tables[(*iter)->GetMaterial()].GetNewE(EejectMod, penetration, dist_traveled);
						}
						else{ std::cout << " ERROR! Doing energy loss on ejectile particle with Z == 0???\n"; }
					}
					else{ 
						if(recoil_part.GetZ() > 0){ // Calculate energy loss for the recoil in the detector
							QDC = ErecoilMod - recoil_tables[(*iter)->GetMaterial()].GetNewE(ErecoilMod, penetration, dist_traveled);
						}
						else{ std::cout << " ERROR! Doing energy loss on recoil particle with Z == 0???\n"; }
					}	
				}
				else{ // Do not do energy loss calculations. The particle leaves all of its energy in the detector.
					dist_traveled = penetration*frand(); // The fraction of the detector which the particle travels through before interacting
					if(proc_eject){ QDC = EejectMod; } // The ejectile may leave any portion of its energy inside the detector
					else{ QDC = ErecoilMod; } // The recoil may leave any portion of its energy inside the detector
				}

				// If particle originates outside of the detector, add the flight path to the first encountered
				// detector face. Otherwise, if the particle originates inside the detector (i.e. a detector at
				// the origin), the distance traveled through the detector is already the total flight path and
				// the time of flight is irrelevant.
				tof = 0.0;
				if(fpath1 >= 0.0 && fpath2 >= 0.0){
					dist_traveled += fpath1;
					
					// Calculate the particle ToF (ns)
					if(proc_eject){ tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMejectMeV()/EejectMod); }
					else{ tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMrecoilMeV()/ErecoilMod); }

					// Smear tof due to PIXIE resolution
					tof += rndgauss0(timeRes); 
				}

				// Get the local coordinates of the intersection point.
				(*iter)->GetLocalCoords(HitDetect1, hit_x, hit_y, hit_z);
			
				// Calculate the lab angles of the detector intersection point. Ignore normalization, we're
				// going to throw away R anyway.
				if(proc_eject){ Cart2Sphere(HitDetect1, EjectSphere); }
				else{ Cart2Sphere(HitDetect1, RecoilSphere); }
			
				// Calculate the hit detection point in 3d space. This point will lie along the vector pointing
				// from the origin to the point where the ray intersects a detector and takes finite range of a
				// particle in a material into account.
				HitDetect1.Normalize();
				HitDetect1 = HitDetect1*dist_traveled;
			
				// Main output
				if(proc_eject){
					EJECTdata.Append(HitDetect1.axis[0], HitDetect1.axis[1], HitDetect1.axis[2], EjectSphere.axis[1]*rad2deg,
									 EjectSphere.axis[2]*rad2deg, QDC, tof*(1E9), rdata.Eeject, hit_x, hit_y, hit_z, (*iter)->GetLoc(), false);
								
					// Adjust the ejectile energy to take energy loss into account. 
					EejectMod = EejectMod - QDC;
				}
				else{
					RECOILdata.Append(HitDetect1.axis[0], HitDetect1.axis[1], HitDetect1.axis[2], RecoilSphere.axis[1]*rad2deg,
									  RecoilSphere.axis[2]*rad2deg, QDC, tof*(1E9), rdata.Erecoil, hit_x, hit_y, hit_z, (*iter)->GetLoc(), false);
				
					// Adjust the recoil energy to take energy loss into account. 
					ErecoilMod = ErecoilMod - QDC;
				}
			} // if(hit)
			
			// Check if we need to process the recoil in this detector.
			if(proc_eject){ 
				proc_eject = false; 
				if((*iter)->IsRecoilDet()){ goto process; }
			}
		} // for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++)
		if(InCoincidence){ // We require coincidence between ejectiles and recoils 
			if(EJECTdata.eject_mult > 0 && RECOILdata.recoil_mult > 0){ 
				if(!flag){ flag = true; }
				if(WriteReaction){ // Set some extra reaction data variables.
					REACTIONdata.Append(rdata.Ereact, rdata.Eeject, rdata.Erecoil, rdata.comAngle*rad2deg, rdata.state,
					                    lab_beam_interaction.axis[0], lab_beam_interaction.axis[1], lab_beam_interaction.axis[2],
					                    lab_beam_stragtraject.axis[0], lab_beam_stragtraject.axis[1], lab_beam_stragtraject.axis[2]);
				}
				if(bgPerDetection){ backgroundWait = backgroundRate; }
				VIKARtree->Fill(); 
				Ndetected++;
			}
		}
		else{ // Coincidence is not required between reaction particles
			if(EJECTdata.eject_mult > 0 || RECOILdata.recoil_mult > 0){ 
				if(!flag){ flag = true; }
				if(WriteReaction){
					REACTIONdata.Append(rdata.Ereact, rdata.Eeject, rdata.Erecoil, rdata.comAngle*rad2deg, rdata.state,
					                    lab_beam_interaction.axis[0], lab_beam_interaction.axis[1], lab_beam_interaction.axis[2],
					                    lab_beam_stragtraject.axis[0], lab_beam_stragtraject.axis[1], lab_beam_stragtraject.axis[2]);
				}
				if(bgPerDetection){ backgroundWait = backgroundRate; }
				VIKARtree->Fill(); 
				Ndetected++;
			}
		}
		
		// Zero all output data structures.
		EJECTdata.Zero();
		RECOILdata.Zero();
		if(WriteReaction){ REACTIONdata.Zero(); }
	} // Main simulation loop
	// ==  ==  ==  ==  ==  ==  == 
	
	// Information output and cleanup
	std::cout << "\n ------------- Simulation Complete --------------\n";
	std::cout << " Simulation Time: " << (float)(clock()-timer)/CLOCKS_PER_SEC << " seconds\n"; 
	std::cout << " Total MC Events: " << Nreactions << "\n";
	std::cout << " Total Detector Hits: " << NdetHit << " (" << NdetHit*100.0/Nreactions << "%)\n";
	if(beam_stopped > 0 || eject_stopped > 0 || recoil_stopped > 0){
		std::cout << " Particles Stopped in Target:\n";
		if(beam_stopped > 0){ std::cout << "  Beam: " << beam_stopped << " (" << 100.0*beam_stopped/Nsimulated << "%)\n"; }
		if(eject_stopped > 0){ std::cout << "  Ejectiles: " << eject_stopped << " (" << 100.0*eject_stopped/Nsimulated << "%)\n"; }
		if(recoil_stopped > 0){ std::cout << "  Recoils: " << recoil_stopped << " (" << 100.0*recoil_stopped/Nsimulated << "%)\n"; }
	}
	if(SupplyRates){ std::cout << " Beam Time: " << Nsimulated/BeamRate << " seconds\n"; }
	
	file->cd();
	VIKARtree->Write();
	if(DEBUGtree){ DEBUGtree->Write(); }
	
	std::cout << "  Wrote file " << output_fname << "\n";
	std::cout << "   Wrote " << VIKARtree->GetEntries() << " tree entries for VIKAR\n";
	if(WriteDebug){ std::cout << "   Wrote " << DEBUGtree->GetEntries() << " tree entries for DEBUG\n"; }
	file->Close();
	
	delete file;
	delete[] materials;
	if(eject_part.GetZ() > 0){ delete[] eject_tables; }
	if(recoil_part.GetZ() > 0){ delete[] recoil_tables; }
	delete[] ExRecoilStates;
	delete[] totXsect;

	for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
		delete *iter;
	}
	vandle_bars.clear();
	
	return 0;
} 
