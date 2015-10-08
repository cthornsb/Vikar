// vikar.h
// Converted by FortranConvert v0.1
// Wed Feb 12 19:55:17 2014

#include <fstream>
#include <iostream>
#include <time.h>

#include "TFile.h"
#include "TTree.h"

#include "vikar_core.h"
#include "kindeux.h"
#include "materials.h"
#include "detectors.h"
#include "Structures.h"

#define VERSION "1.21b"

struct debugData{
	double var1, var2, var3;
	
	void Set(double v1, double v2, double v3){
		var1 = v1; var2 = v2; var3 = v3;
	}
};

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
	Planar *vandle_bars = NULL; // Array of Planar detectors

	RangeTable beam_targ; // Range table for beam in target
	//RangeTable *eject_targ = NULL; // Pointer to the range table for ejectile in target
	//RangeTable *recoil_targ = NULL; // Pointer to the range table for recoil in target
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
	std::string output_fname_prefix = "VIKAR";
	
	// Input/output variables
	int face1, face2; // Detect which detector faces are hit
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
	bool InCoincidence = true;
	bool WriteReaction = false;
	bool WriteDebug = false;
	bool PerfectDet = true;
	bool ADists = false;
	bool SupplyRates = false;
	bool BeamFocus = false;
	bool DoRutherford = false;
	bool CylindricalBeam = true;

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
		if(argc >= 3){ output_fname_prefix = std::string(argv[2]); }
		
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
				targ.SetZ((double)atof(input.c_str()));
				std::cout << "  Target-Z: " << targ.GetZ() << std::endl;
			}
			else if(count == 4){ 
				targ.SetA((double)atof(input.c_str()));
				std::cout << "  Target-A: " << targ.GetA() << std::endl;
			}
			else if(count == 5){ 
				eject_part.SetZ(atof(input.c_str()));
				recoil_part.SetZ(beam_part.GetZ() + targ.GetZ() - eject_part.GetZ());
				std::cout << "  Ejectile-Z: " << eject_part.GetZ() << std::endl;
			}
			else if(count == 6){ 
				eject_part.SetA(atof(input.c_str())); 
				recoil_part.SetA(beam_part.GetA() + targ.GetA() - eject_part.GetA()); 
				std::cout << "  Ejectile-A: " << eject_part.GetA() << std::endl;
				std::cout << "  Recoil-Z: " << recoil_part.GetZ() << std::endl;
				std::cout << "  Recoil-A: " << recoil_part.GetA() << std::endl;
			}
			else if(count == 7){ 
				Ebeam0 = atof(input.c_str());  
				std::cout << "  Beam Energy: " << Ebeam0 << " MeV\n";
			}
			else if(count == 8){ 
				SetBool(input, "  Cylindrical Beam", CylindricalBeam);
				getline(input_file, input); input = Parse(input);
				beamspot = atof(input.c_str());
				if(CylindricalBeam){ std::cout << "  Beam Spot Diameter: " << beamspot << " mm\n"; }
				else{ std::cout << "  Beam Spot FWHM: " << beamspot << " mm\n"; }
				beamspot = beamspot/1000.0; // in meters
			}
			else if(count == 9){ 
				beamAngdiv = atof(input.c_str());  
				std::cout << "  Beam Angular Divergence: " << beamAngdiv << " degrees\n";
				beamAngdiv *= deg2rad; // in radians
			}
			else if(count == 10){ 
				beamEspread = atof(input.c_str());  
				std::cout << "  Beam Energy Spread: " << beamEspread << " MeV\n";
			}
			else if(count == 11){ 
				gsQvalue = atof(input.c_str());  
				std::cout << "  G.S. Q-Value: " << gsQvalue << " MeV\n";
			}
			else if(count == 12){ 
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
			else if(count == 13){ 
				// Angular distribution information
				if(SetBool(input, "  Supply Angular Distributions", ADists)){
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
					
					// Supply beam rate information
					getline(input_file, input); input = Parse(input);
					if(SetBool(input, "  Calculate Rates", SupplyRates)){ 
						getline(input_file, input); input = Parse(input);
						BeamRate = atof(input.c_str());
						std::cout << "   Beam Rate: " << BeamRate << " pps\n";
					}
					else{ BeamRate = 0.0; }
				}
				else{ SupplyRates = false; }
			}
			else if(count == 14){
				// Target material
				targ_mat_name = input;
			}
			else if(count == 15){ 
				// Target thickness
				targ.SetThickness((double)atof(input.c_str()));
				std::cout << "  Target Thickness: " << targ.GetThickness() << " mg/cm^2\n";			
			}
			else if(count == 16){ 
				// Target angle wrt beam axis
				targ.SetAngle((double)atof(input.c_str())*deg2rad);
				std::cout << "  Target Angle: " << targ.GetAngle()*rad2deg << " degrees\n";
			}
			else if(count == 17){ 
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
			else if(count == 18){ 
				// Load detector setup from a file
				det_fname = "./detectors/" + input;
				std::cout << "  Detector Setup Filename: " << det_fname << std::endl;
			}
			else if(count == 19){ 
				// Desired number of detections
				Nwanted = atol(input.c_str());
				std::cout << "  Desired Detections: " << Nwanted << std::endl; 
			}
			else if(count == 20){
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
			else if(count == 21){
				// Require ejectile and recoil particle coincidence?
				SetBool(input, "  Require Particle Coincidence", InCoincidence);
			}
			else if(count == 22){
				// Write Reaction data to file?
				SetBool(input, "  Write Reaction Info", WriteReaction);
			}
			else if(count == 23){
				// Write Debug data to file?
				SetBool(input, "  Write Debug Info", WriteDebug);
			}
			
			count++;
		}
		
		input_file.close();
		if(count <= 23){ std::cout << " Warning! The input file is invalid. Check to make sure input is correct\n"; }
	}
	else{
		std::cout << " Error! Missing required variable\n";
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
	for(int i = 0; i < Ndet; i++){
		if(!IsInVector(vandle_bars[i].GetMaterialName(), needed_materials)){
			needed_materials.push_back(vandle_bars[i].GetMaterialName());
		}
		
		if(vandle_bars[i].IsEjectileDet()){ NdetEject++; }
		if(vandle_bars[i].IsRecoilDet()){ NdetRecoil++; }
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
		std::cout << "  Target Radiation Length: " << targ.GetRadLength() << " mg/cm^2\n";
	
		// Calculate the stopping power table for the beam particles in the target
		if(beam_part.GetZ() > 0){ // The beam is a charged particle (not a neutron)
			std::cout << "\n Calculating range table for beam in " << materials[targ_mat_id].GetName() << "...";
			beam_targ.Init(100, 0.1, (Ebeam0+2*beamEspread), beam_part.GetZ(), beam_part.GetA()*amu2mev, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
	}

	// Calculate the stopping power table for the ejectiles in the materials
	if(eject_part.GetZ() > 0){ // The ejectile is a charged particle (not a neutron)
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating ejectile range table for " << materials[i].GetName() << "...";
			eject_tables[i].Init(100, 0.1, (Ebeam0+2*beamEspread), eject_part.GetZ(), eject_part.GetA()*amu2mev, &materials[i]);
			std::cout << " Done!\n";
		}
		//eject_targ = &eject_tables[targ_mat_id]; // Table for ejectile in target
	}


	// Calculate the stopping power table for the recoils in the materials
	if(recoil_part.GetZ() > 0){ // The recoil is a charged particle (not a neutron)
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating recoil range table for " << materials[i].GetName() << "...";
			recoil_tables[i].Init(100, 0.1, (Ebeam0+2*beamEspread), recoil_part.GetZ(), recoil_part.GetA()*amu2mev, &materials[i]);
			std::cout << " Done!\n";
		}
		//recoil_targ = &recoil_tables[targ_mat_id]; // Table for recoil in target
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
	for(int i = 0; i < Ndet; i++){ // Set the detector material for energy loss calculations
		for(unsigned int j = 0; j < num_materials; j++){
			if(vandle_bars[i].GetMaterialName() == materials[j].GetName()){
				if((vandle_bars[i].IsRecoilDet() && recoil_part.GetZ() > 0) || (vandle_bars[i].IsEjectileDet() && eject_part.GetZ() > 0)){ 
					// Only set detector to use material if the particle it is responsible for detecting has
					// a Z greater than zero. Particles with Z == 0 will not have calculated range tables
					// and thus cannot use energy loss considerations.
					vandle_bars[i].SetMaterial(j);  
				}
				break; 
			}
		}
	}

	//std::cout << "\n ==  ==  ==  ==  == \n\n";

	if(ADists && !DoRutherford){ 
		std::cout << "\n Loading state angular distribution files...\n";
		if(kind.SetDist(AngDist_fname, materials[targ_mat_id].GetTotalElements(), materials[targ_mat_id].GetDensity(), BeamRate)){
			// Successfully set the angular distributions
			std::cout << " Successfully loaded angular distributions\n";
			kind.Print();
		}
		else{
			std::cout << "  Warning! Could not properly initialize distributions.\n";
			std::cout << "  Note: Setting all energy states to isotropic!\n";
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
	TFile *file = new TFile("VIKAR.root", "RECREATE");
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
			for(int bar = 0; bar < Ndet; bar++){
				if(!(vandle_bars[bar].IsEjectileDet() || vandle_bars[bar].IsRecoilDet())){ continue; } // This detector cannot detect particles

				// Select the "tof" of the background event
				tof = (double)frand(0, detWindow)*(1E-9);
			
				// Select a random point insde the detector
				vandle_bars[bar].GetRandomPointInside(temp_vector);
				Cart2Sphere(temp_vector_sphere);
			
				// Calculate the apparent energy of the particle using the tof
				if(vandle_bars[bar].IsEjectileDet()){
					double dummyE = 0.5*kind.GetMejectMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					//if(vandle_bars[bar].GetType() == "vandle"){ dummyE = MeV2MeVee(dummyE); }
					EJECTdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], temp_vector_sphere.axis[1]*rad2deg,
									 temp_vector_sphere.axis[2]*rad2deg, dummyE, tof*(1E9), 0.0, 0.0, 0.0, bar, true);
					VIKARtree->Fill(); 
					EJECTdata.Zero();
				}
				else if(vandle_bars[bar].IsRecoilDet()){
					double dummyE = 0.5*kind.GetMrecoilMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					RECOILdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], RecoilSphere.axis[1]*rad2deg,
									  RecoilSphere.axis[2]*rad2deg, dummyE, tof*(1E9), 0.0, 0.0, 0.0, bar, true);
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
						Vector3 dum1, dum2;
						int temp1, temp2;
						double tempx, tempy, tempz;
						targ.GetPlanar()->IntersectPrimitive(Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, 1.0), dum1, dum2, temp1, temp2, tempx, tempy, tempz);
						std::cout << "  Front face intersect = (" << dum1.Dump() << ")\n";
						std::cout << "  Back face intersect = (" << dum2.Dump() << ")\n";
						std::cout << "  Target thickness: " << (dum2-dum1).Length() << " m\n";
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

			EejectMod = rdata.Eeject;
			ErecoilMod = rdata.Erecoil;

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

			// Calculate the energy loss for the ejectile and recoils in the target
			/*if(use_target_eloss){
				if(eject_part.GetZ() > 0){
					range_eject = eject_targ->GetRange(Eeject);
				}
				if(recoil_part.GetZ() > 0){
					range_recoil = recoil_targ->GetRange(Erecoil);
				}
			}*/
		
			if(WriteDebug){ 
				DEBUGdata.Set(Ejectile.axis[0], Ejectile.axis[1], Ejectile.axis[2]);
				DEBUGtree->Fill();
			}
		}

		// Process the reaction products
		for(int bar = 0; bar < Ndet; bar++){
			if(vandle_bars[bar].IsEjectileDet()){ // This detector can detect ejectiles
				if(EejectMod <= 0.0){ continue; } // There still may be recoil detectors left
				proc_eject = true; 
			}
			else if(vandle_bars[bar].IsRecoilDet()){ // This detector can detect recoils
				if(ErecoilMod <= 0.0){ continue; } // There still may be ejectile detectors left
				proc_eject = false; 
			}
			else{ continue; } // This detector cannot detect particles, so skip it

process:			
			if(proc_eject){ hit = vandle_bars[bar].IntersectPrimitive(lab_beam_interaction, Ejectile, HitDetect1, HitDetect2, face1, face2, hit_x, hit_y, hit_z); }
			else{ hit = vandle_bars[bar].IntersectPrimitive(lab_beam_interaction, Recoil, HitDetect1, HitDetect2, face1, face2, hit_x, hit_y, hit_z); }
			
			// If a geometric hit was detected, process the particle
			if(hit){
				NdetHit++; 
				
				// The time of flight is the time it takes the particle to traverse the distance
				// from the target to the intersection point inside the detector
				fpath1 = HitDetect1.Length(); // Distance from reaction to first intersection point
				fpath2 = HitDetect2.Length(); // Distance from reaction to second intersection point
				temp_vector = (HitDetect2-HitDetect1); // The vector pointing from the first intersection point to the second
				penetration = temp_vector.Length(); // Total distance traveled through detector	
				dist_traveled = 0.0;

				if(vandle_bars[bar].UseMaterial()){ // Do energy loss and range considerations
					if(proc_eject){
						if(eject_part.GetZ() > 0){ // Calculate energy loss for the ejectile in the detector
							QDC = EejectMod - eject_tables[vandle_bars[bar].GetMaterial()].GetNewE(EejectMod, penetration);
						}
						else{ std::cout << " ERROR! Doing energy loss on ejectile particle with Z == 0???\n"; }
					}
					else{ 
						if(recoil_part.GetZ() > 0){ // Calculate energy loss for the recoil in the detector
							QDC = ErecoilMod - recoil_tables[vandle_bars[bar].GetMaterial()].GetNewE(ErecoilMod, penetration);
						}
						else{ std::cout << " ERROR! Doing energy loss on recoil particle with Z == 0???\n"; }
					}	
				}
				else{ // Do not do energy loss calculations
					dist_traveled = penetration*frand(); // The fraction of the detector which the particle travels through before interacting
					if(proc_eject){ QDC = EejectMod*frand(); } // The ejectile may leave any portion of its energy inside the detector
					else{ QDC = ErecoilMod*frand(); } // The recoil may leave any portion of its energy inside the detector
				}

				// Calculate the total distance traveled and the interaction point inside the detector
				// Ensure that we use the intersection point on the side facing the target
				if(fpath1 <= fpath2){ 
					dist_traveled += fpath1; 
					temp_vector = lab_beam_interaction + HitDetect1 + temp_vector*penetration;
				}
				else{ 
					dist_traveled += fpath2; 
					temp_vector = lab_beam_interaction + HitDetect2 - temp_vector*penetration;
				}
			
				// Calculate the particle ToF (ns)
				tof = 0.0;
				if(proc_eject){ tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMejectMeV()/EejectMod); }
				else{ tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMrecoilMeV()/ErecoilMod); }

				// Smear tof due to PIXIE resolution
				tof += rndgauss0(timeRes); 
		
				// Main output
				// X(m) Y(m) Z(m) LabTheta(deg) LabPhi(deg) QDC(MeV) ToF(ns) Bar# Face# HitX(m) HitY(m) HitZ(m)
				if(proc_eject){
					double dummyE = 0.5*kind.GetMejectMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					if(vandle_bars[bar].GetType() == "vandle"){ dummyE = MeV2MeVee(dummyE); }
					Cart2Sphere(temp_vector, EjectSphere); // Ignore normalization, we're going to throw away R anyway
					EJECTdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], EjectSphere.axis[1]*rad2deg,
									 EjectSphere.axis[2]*rad2deg, dummyE, tof*(1E9), hit_x, hit_y, hit_z, bar, false);
				}
				else{
					double dummyE = 0.5*kind.GetMrecoilMeV()*dist_traveled*dist_traveled/(c*c*tof*tof);
					Cart2Sphere(temp_vector, RecoilSphere); // Ignore normalization, we're going to throw away R anyway
					RECOILdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], RecoilSphere.axis[1]*rad2deg,
									  RecoilSphere.axis[2]*rad2deg, dummyE, tof*(1E9), hit_x, hit_y, hit_z, bar, false);
				}
				
				// Adjust the particle energies to take energy loss into account
				if(proc_eject){ EejectMod = EejectMod - QDC; }
				else{ ErecoilMod = ErecoilMod - QDC; }
			} // if(hit)
			
			if(proc_eject){ 
				proc_eject = false; 
				if(vandle_bars[bar].IsRecoilDet()){ goto process; } // Still need to process the recoil in this detector
			}
		} // for(unsigned int bar = 0; bar < Ndet; bar++)
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
	
	std::cout << "  Wrote file " << output_fname_prefix << ".root\n";
	std::cout << "   Wrote " << VIKARtree->GetEntries() << " tree entries for VIKAR\n";
	if(WriteDebug){ std::cout << "   Wrote " << DEBUGtree->GetEntries() << " tree entries for DEBUG\n"; }
	file->Close();
	
	delete file;
	delete[] materials;
	if(eject_part.GetZ() > 0){ delete[] eject_tables; }
	if(recoil_part.GetZ() > 0){ delete[] recoil_tables; }
	delete[] vandle_bars;
	delete[] ExRecoilStates;
	delete[] totXsect;
	
	return 0;
} 
