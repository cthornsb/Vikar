// vikar.h
// Converted by FortranConvert v0.1
// Wed Feb 12 19:55:17 2014

#include <fstream>
#include <iostream>
#include <time.h>

#include "../include/vikar_core.h"
#include "../include/planar.h"

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
	Kindeux kind;
	Efficiency bar_eff;
	
	// temporary variables
	double hit_x, hit_y, hit_z;
	
	Vector3 Ejectile, Gamma, HitDetect, HitDetectSphere, RecoilSphere;
	
	// Energy loss varibales
	double tof;
	unsigned short NtargElements; 
	unsigned short **targComp; 
	double tgtdens;
	double thickness, Zdepth, targ_depth, targ_angle; 
	double Ereact;
	double AveTargZ, AveTargA, totalTarg; 
	double targRadL;
	
	// Physics Variables
	unsigned short NRecoilStates;
	std::vector<std::string> AngDist_fname; 
	double *ExRecoilStates, gsQvalue; 
	double *totXsect; 
	double Ebeam, Ebeam0, Erecoil, Eeject, Abeam, Zbeam, Atarg;
	double Arecoil, Zrecoil, Aeject, Zeject;
	double Ztarg;
	
	// Beam variables
	double thetaBeam, phiBeam; // Incident beam particle direction
	double thetaReact, phiReact; // Direction of beam particle at reaction (after straggling)
	double beamspot, beamEspread; // Beamspot size FWHM (mm)

	double timeRes = 3E-9;
	double energyRes;
	double BeamRate;
	
	// Detector variables
	std::string det_fname; 
	Planar *vandle_bars;
	
	std::string output_fname_prefix = "VIKAR";
	
	// Input/output variables
	short face;
	int Ndetected, Nwanted;
	int Nsimulated, NbarHit;
	int Nreactions;
	clock_t timer;

	bool SimGamma, WriteRecoil, WriteDebug;

	unsigned short Ndet;
	bool PerfectDet, DetSetup, ADists, SupplyRates;
	short idummy;

	// Error reporting variables
	float version;
	std::string SRIM_fName, SRIM_fName_beam, SRIM_fName_targ_eject, SRIM_fName_targ_recoil;
	std::string SRIM_fName_det_eject, SRIM_fName_det_recoil;

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

	std::cout << "\n VANDLE VIKAR 1.0\n"; 
	std::cout << " ==  ==  ==  ==  == \n\n"; 
	
	std::cout << " Welcome to NewVIKAR, the Virtual Instrumentation for Kinematics\n"; 
	std::cout << " And Reactions program, optimized for use with VANDLE bars\n\n";

	std::cout << " How about a nice cup of tea?\n"; 

	std::ifstream input_file;
	if(argc >= 2){
		// Read an input file
		input_file.open(argv[1]);
		std::cout << " No? Well let me just load the input file then\n\n ==  ==  ==  ==  == \n";
		if(!input_file.good()){
			std::cout << " Error: Problem loading the input file\n";
			return 1;
		}
		
		// Specify the name of the output files
		if(argc >= 3){ output_fname_prefix = std::string(argv[2]); }
		
		unsigned short count = 0;
		std::string input;
		std::cout << "\n Reading from file " << argv[1] << std::endl;
		while(true){
			input_file >> input;
			if(input_file.eof()){ break; }
			if(count == 0){ 
				version = atof(input.c_str()); 
				std::cout << "  Version: " << version << std::endl;
				if(version != 1.0){ std::cout << "   Warning! This input file has the wrong version number. Check to make sure input is correct\n"; }
			}
			else if(count == 1){ 
				Zbeam = atof(input.c_str());  
				std::cout << "  Beam-Z: " << Zbeam << std::endl;
			}
			else if(count == 2){ 
				Abeam = atof(input.c_str());  
				std::cout << "  Beam-A: " << Abeam << std::endl;
			}
			else if(count == 3){ 
				Ztarg = atof(input.c_str());  
				std::cout << "  Target-Z: " << Ztarg << std::endl;
			}
			else if(count == 4){ 
				Atarg = atof(input.c_str());  
				std::cout << "  Target-A: " << Atarg << std::endl;
			}
			else if(count == 5){ 
				Zeject = atof(input.c_str());
				Zrecoil = Zbeam + Ztarg - Zeject;  
				std::cout << "  Ejectile-Z: " << Zeject << std::endl;
			}
			else if(count == 6){ 
				Aeject = atof(input.c_str()); 
				Arecoil = Abeam + Atarg - Aeject; 
				std::cout << "  Ejectile-A: " << Aeject << std::endl;
				std::cout << "  Recoil-Z: " << Zrecoil << std::endl;
				std::cout << "  Recoil-A: " << Arecoil << std::endl;
			}
			else if(count == 7){ 
				Ebeam0 = atof(input.c_str());  
				std::cout << "  Beam Energy: " << Ebeam0 << " MeV\n";
			}
			else if(count == 8){ 
				beamspot = atof(input.c_str());  
				std::cout << "  Beam Spot Size: " << beamspot << " mm\n";
				beamspot *= 0.1;
			}
			else if(count == 9){ 
				beamEspread = atof(input.c_str());  
				std::cout << "  Beam Spread: " << beamEspread << " MeV\n";
			}
			else if(count == 10){ 
				gsQvalue = atof(input.c_str());  
				std::cout << "  G.S. Q-Value: " << gsQvalue << " MeV\n";
			}
			else if(count == 11){ 
				// Recoil excited state information
				NRecoilStates = atos(input.c_str()) + 1;
				std::cout << "  No. Excited States: " << NRecoilStates-1 << std::endl;
				ExRecoilStates = new double[NRecoilStates]; ExRecoilStates[0] = 0.0;
				totXsect = new double[NRecoilStates]; totXsect[0] = 0.0;
				std::cout << "   Recoil Ground State: 0.0 MeV\n";
				for(unsigned short i = 1; i < NRecoilStates; i++){
					input_file >> input;
					ExRecoilStates[i] = atof(input.c_str());
					std::cout << "   Recoil State " << i+1 << ": " << ExRecoilStates[i] << " MeV\n";
					totXsect[i] = 0.0;
				}
			}
			else if(count == 12){ 
				// Angular distribution information
				idummy = atos(input.c_str());
				if(idummy == 1){ ADists = true; }
				else{ ADists = false; }
				std::cout << "  Supply Angular Distributions: ";
				if(ADists){
					std::cout << "Yes\n";
					for(unsigned short i = 0; i < NRecoilStates; i++){
						input_file >> input;
						AngDist_fname.push_back(input);
						if(i == 0){ std::cout << "   Distribution for ground state: " << AngDist_fname[i] << std::endl; }
						else{ std::cout << "   Distribution for state " << i+1 << ": " << AngDist_fname[i] << std::endl; }
					}
					
					// Supply beam rate information
					input_file >> input;
					idummy = atos(input.c_str());
					if(idummy == 1){ SupplyRates = true; }
					else{ SupplyRates = false; }
					std::cout << "  Calculate Rates: ";
					if(SupplyRates){ 
						std::cout << "Yes\n"; 
						input_file >> input;
						BeamRate = atof(input.c_str());
						std::cout << "   Beam Rate: " << BeamRate << " pps\n";
					}
					else{ 
						std::cout << "No\n"; 
						BeamRate = 0.0;
					}
				}
				else{ 
					std::cout << "No\n"; 
					SupplyRates = false;
				}
			}
			else if(count == 13){ 
				// Target information
				NtargElements = atos(input.c_str());
				std::cout << "  No. Target Elements: " << NtargElements << std::endl;
				targComp = new unsigned short *[NtargElements];
				for(unsigned short i = 0; i < NtargElements; i++){
					targComp[i] = new unsigned short[3];
					input_file >> input; targComp[i][0] = atos(input.c_str()); 
					input_file >> input; targComp[i][1] = atos(input.c_str());  
					input_file >> input; targComp[i][2] = atos(input.c_str()); 
					std::cout << "   Element " << i+1 << ": " << targComp[i][2] << " per molecule of Z = " << targComp[i][0] << ", A = " << targComp[i][1] << std::endl;
				}
				
				// Calculate the radiation length for the target material
				// for use in angular straggling calculations later.
				// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
				AveTargZ = 0.0; 
				AveTargA = 0.0; 
				totalTarg = 0.0; 
				targRadL = 0.0; 
				for(unsigned short i = 0; i < NtargElements; i++){
					AveTargZ += targComp[i][0]*targComp[i][2]; 
					AveTargA += targComp[i][1]*targComp[i][2]; 
					totalTarg += targComp[i][2]; 
					targRadL += ((targComp[i][1]*targComp[i][2]/AveTargA)/radlength(targComp[i][1],targComp[i][0])); 
				} 
	
				targRadL = 1.0/targRadL; 
				AveTargZ = AveTargZ/totalTarg; 
				AveTargA = AveTargA/totalTarg;
			}
			else if(count == 14){ 
				// Target thickness
				thickness = atof(input.c_str());  
				std::cout << "  Target Thickness: " << thickness << " mg/cm^2\n";
			}
			else if(count == 15){ 
				// Target density
				tgtdens = atof(input.c_str());  
				std::cout << "  Target Density: " << tgtdens << " g/cm^3\n";
			}
			else if(count == 16){ 
				// Target angle wrt beam axis
				targ_angle = atof(input.c_str());  
				std::cout << "  Target Angle: " << targ_angle << " degrees\n";
				targ_angle *= deg2rad;
			}
			else if(count == 17){ 
				// Load the small, medium, and large bar efficiencies
				// Efficiency index 0 is the underflow efficiency (for energies below E[0])
				// Efficiency index N is the overflow efficiency (for energies greater than E[N])
				idummy = atos(input.c_str());
				if(idummy == 1){ PerfectDet = true; }
				else{ PerfectDet = false; }
				std::cout << "  Perfect Detector: ";
				if(!PerfectDet){ 
					// Load small bar efficiency data
					std::cout << "No\n";
					input_file >> input; 
					std::cout << "   Found " << bar_eff.ReadSmall(input.c_str()) << " small bar data points in file " << input << "\n";
					
					// Load medium bar efficiency data
					input_file >> input; 
					std::cout << "   Found " << bar_eff.ReadMedium(input.c_str()) << " medium bar data points in file " << input << "\n";
					
					// Load large bar efficiency data
					input_file >> input; 
					std::cout << "   Found " << bar_eff.ReadLarge(input.c_str()) << " large bar data points in file " << input << "\n";
				}
				else{ std::cout << "Yes\n"; }
			}
			else if(count == 18){ 
				// Supply detector setup file?
				idummy = atos(input.c_str());
				if(idummy == 1){ DetSetup = true; }
				else{ DetSetup = false; }
				std::cout << "  Detector Setup File: ";
				if(DetSetup){
					// Load detector setup from a file
					std::cout << "Yes\n";
					input_file >> det_fname;
					std::cout << "   Path: " << det_fname << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 19){ 
				// Desired number of detections
				Nwanted = atol(input.c_str());
				std::cout << "  Desired Detections: " << Nwanted << std::endl; 
			}
			else if(count == 20){
				// Simulate prompt gamma flash?
				idummy = atos(input.c_str());
				if(idummy == 1){ SimGamma = true; }
				else{ SimGamma = false; }
				std::cout << "  Detect Prompt Gammas: ";
				if(SimGamma){ std::cout << "Yes\n"; }
				else{ std::cout << "No\n"; }
			}
			else if(count == 21){
				// Write Recoil data to file?
				idummy = atos(input.c_str());
				if(idummy == 1){ WriteRecoil = true; }
				else{ WriteRecoil = false; }
				std::cout << "  Write Recoil: ";
				if(WriteRecoil){ std::cout << "Yes\n"; }
				else{ std::cout << "No\n"; }
			}
			else if(count == 22){
				// Write Debug data to file?
				idummy = atos(input.c_str());
				if(idummy == 1){ WriteDebug = true; }
				else{ WriteDebug = false; }
				std::cout << "  Write Debug Info: ";
				if(WriteDebug){ std::cout << "Yes\n"; }
				else{ std::cout << "No\n"; }
			}
			
			count++;
		}
		
		input_file.close();
		if(count < 22){
			std::cout << " Error: The input file is invalid\n";
			return 1;
		}
	}
	else{
		std::cout << " Error: Missing required variable\n";
		return 1;
	}
	
	std::cout << "\n ==  ==  ==  ==  == \n\n";

	// Read VIKAR detector setup file or manually setup simple systems
	if(DetSetup){
		std::cout << " Reading in NewVIKAR detector setup file...\n";

		std::ifstream detfile(det_fname.c_str());
		if(!detfile.good()){
			std::cout << " Warning: Failed to load detector file " << det_fname << std::endl;
			return 1;
		}

		std::vector<NewVIKARDet> detectors;
		std::string bar_type;
		float values[9];

		while(true){
			for(unsigned short i = 0; i < 6; i++){ 
				detfile >> values[i];
			}

			// Set the size of the bar
			detfile >> bar_type;
			if(!(bar_type == "small" || bar_type == "medium" || bar_type == "large")){
				for(unsigned short i = 6; i < 9; i++){ detfile >> values[i]; }
			}

			detectors.push_back(NewVIKARDet(values, bar_type));
			if(detfile.eof()){ break; }
		}	
		detfile.close();

		// Generate the Planar bar array
		vandle_bars = new Planar[detectors.size()];
		for(unsigned short i = 0; i < detectors.size(); i++){
			// Set the size of the bar
			if(detectors[i].bar_type == "small"){ vandle_bars[i].SetSmall(); }
			else if(detectors[i].bar_type == "medium"){ vandle_bars[i].SetMedium(); }
			else if(detectors[i].bar_type == "large"){ vandle_bars[i].SetLarge(); }
			else{ vandle_bars[i].SetSize(detectors[i].data[6],detectors[i].data[7],detectors[i].data[8]); }

			vandle_bars[i].SetPosition(detectors[i].data[0],detectors[i].data[1],detectors[i].data[2]); // Set the x,y,z position of the bar
			vandle_bars[i].SetBarRotation(detectors[i].data[3],detectors[i].data[4],detectors[i].data[5]); // Set the 3d rotation of the bar
		}

		Ndet = detectors.size();

		// Report on how many detectors were read in
		std::cout << " Found " << Ndet << " VANDLE detectors in file " << det_fname << std::endl;

		// Check there's at least 1 detector!
		if(Ndet < 1){
			std::cout << " Error: Found no detectors. Check that the filename is correct\n"; 
			return 1;
		}
	}
	else{
		// Manual detector setup (fix later)
		std::cout << " Not Implemented!\n";
		return 1;
	}

	//std::cout << "\n ==  ==  ==  ==  == \n\n";
	std::cout << "\n Initializing main simulation Kindeux object...\n";

	// Initialize kinematics object
	kind.Initialize(Abeam, Atarg, Arecoil, Aeject, gsQvalue, NRecoilStates, ExRecoilStates, tgtdens);
	if(ADists){ 
		std::cout << " Loading state angular distribution files...\n";
		if(kind.SetDist(AngDist_fname, totalTarg, BeamRate)){
			// Successfully set the angular distributions
			kind.Print();
		}
		else{
			std::cout << "  Warning! Could not properly initialize distributions.\n";
			std::cout << "  Note: Setting all energy states to isotropic!\n";
		}
	}

	std::cout << "\n ==  ==  ==  ==  == \n\n";

	//---------------------------------------------------------------------------
	// End of Input Section
	//---------------------------------------------------------------------------
	Ndetected = 0; // Total number of particles detected in VANDLE
	NbarHit = 0; // Total number of particles which collided with a bar
	Nsimulated = 0; // Total number of simulated particles
	Nreactions = 0; // Total number of particles which react with the target

	// Open the output file
	std::ofstream VIKARout((output_fname_prefix+"_main.dat").c_str());
	std::ofstream RecoilOut, DebugOut;
	if(WriteRecoil){ RecoilOut.open((output_fname_prefix+"_recoil.dat").c_str()); }
	if(WriteDebug){ DebugOut.open((output_fname_prefix+"_debug.dat").c_str()); }

	// Begin the simulation
	std::cout << " ---------- Simulation Setup Complete -----------\n"; 
	std::cout << "\n Beginning simulating a bunch of events....\n"; 

	//---------------------------------------------------------------------------
	// The Event Loop
	// ==  ==  ==  ==  ==  ==  == 
	// (Just to make it obvious)
	//---------------------------------------------------------------------------

	int chunk = Nwanted/10;
	float totTime = 0.0;
	char counter = 1;
	bool flag = false;
	bool hit;
	timer = clock();

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
				std::cout << "  Geometric Efficiency: " << NbarHit*100.0/Nreactions << "%\n";
				std::cout << "  Detection Efficiency: " << Ndetected*100.0/Nreactions << "%\n"; 
			}
			if(SupplyRates){ std::cout << "  Beam Time: " << Nsimulated/BeamRate << " seconds\n"; }
			
			std::cout << "  Simulation Time: " << totTime << " seconds\n";
			std::cout << "  Time reamining: " << (totTime/counter)*(10-counter) << " seconds\n";
			counter++; 
		}
		
		Nsimulated++; 

		// Recoil Hit????
		// Set the beam direction to (0,0,1) i.e. along the z axis
		// Set the incident direction for the beam
		thetaBeam = rndgauss0(0.1*deg2rad); 
		phiBeam = frand()*2.0*pi; 

		// Calculate the depth in the target for the reaction
		targ_depth = thickness*frand(); 

		// Calculate thickness traversed by the beam particle (Zdepth),
		// dependent on the target angle and beam particle direction.
		// (use (thickness-targ_depth) as subroutine written for recoil/ejectile)
		targ_thick(thetaBeam,phiBeam,thickness,(thickness-targ_depth),targ_angle,Zdepth); 

		// Calculate the beam particle energy, varied with energy spread
		Ebeam = Ebeam0 + rndgauss0(beamEspread); 

		// Calculate the beam particle range in the target
		//range_beam = linear(beam_E[short(Ebeam)],beam_E[short(Ebeam+1.0)],beam_Erange[short(Ebeam)],beam_Erange[short(Ebeam+1.0)],Ebeam); 
		//range_beam = linear(beam_E[short(Ebeam)-1],beam_E[short(Ebeam)],beam_Erange[short(Ebeam)-1],beam_Erange[short(Ebeam)],Ebeam); 

		// Code for calculating energy loss of particle...
		// Find the ranges for energies straddling the beam energy
		//i = 1; 
		//while(beam_Erange[i] < (range_beam-Zdepth)){ i = i + 1; }
		/*for(unsigned short i = 1; i < num_beam_E; i++){
			if(beam_Erange[i] < (range_beam-Zdepth)){ break; }
		}*/

		// Calculate the new energy
		//Ereact = linear(beam_Erange[i-1],beam_Erange[i],beam_E[i-1],beam_E[i],(range_beam-Zdepth)); // FIX!!!!
		Ereact = Ebeam; // Remove beam particle energy loss effects

		// Determine the angle of the beam particle's trajectory at the
		// interaction point, due to angular straggling and the incident trajectory.
		strag_targ(Abeam,Zbeam,Zdepth,thetaBeam,phiBeam,Ebeam,thetaReact,phiReact,targRadL); 

		// the 2 body kinematics routine to generate the ejectile and recoil
		if(WriteRecoil || SimGamma){ 
			// Need to calculate parameters for the recoil
			if(kind.FillVars(Ereact, thetaReact, phiReact, Eeject, Erecoil, HitDetectSphere, RecoilSphere)){ Nreactions++; }
			else{ continue; } // A reaction did not occur
		}
		else{ 
			// Only interested in the ejectile	
			if(kind.FillVars(Ereact, thetaReact, phiReact, Eeject, HitDetectSphere)){ Nreactions++; }
			else{ continue; } // A reaction did not occur
		}
		
		Sphere2Cart(HitDetectSphere, Ejectile); // HitDetectSphere is a unit vector (no need to normalize)
		if(WriteDebug){ DebugOut << HitDetectSphere.axis[1] << "\t" << HitDetectSphere.axis[2] << Eeject << "\n"; }

		// Gammas are indistinguishable from the ejectile (to a VANDLE bar)
		/*if(SimGamma && RecoilState > 0){
			// Simulate gamma decays for excited recoils
			double gamma_theta, gamma_phi;
			UnitRandom(gamma_theta, gamma_phi);
		}*/

		// Process the ejectile
		for(unsigned short bar = 0; bar < Ndet; bar++){
			face = -1;
			face = vandle_bars[bar].FaceIntersect(Ejectile, HitDetect, hit_x, hit_y, hit_z);
			if(face != -1){
				// Geometric hit detected
				NbarHit++;
				hit = true;
				
				// Check for a "true" hit
				if(!PerfectDet){
					// Detector is not perfect, hit is based on the efficiency
					if(vandle_bars[bar].IsSmall()){
						// Use the small bar efficiency data
						if(frand() > bar_eff.GetSmallEfficiency(Eeject)){ hit = false; }
					}
					else if(vandle_bars[bar].IsMedium()){
						// Use the medium bar efficiency data
						if(frand() > bar_eff.GetMediumEfficiency(Eeject)){ hit = false; }
					}
					else if(vandle_bars[bar].IsLarge()){
						// Use the large bar efficiency data
						if(frand() > bar_eff.GetLargeEfficiency(Eeject)){ hit = false; }
					}
				}
				
				if(hit){
					// The neutron hit a bar and was detected
					// The time of flight is the time it takes the neutron to traverse the distance
					// from the target (origin) to the intersection point ON THE FACE of the bar
					tof = HitDetect.Length()*std::sqrt(kind.GetMeject()/(2*Eeject*1.60217657E-13*6.02214129E26)); // Neutron ToF (s)

					// Smear ToF and Energy if the detector is not perfect
					if(!PerfectDet){
						tof += rndgauss0(timeRes); // Smear tof due to VANDLE resolution
						energyRes = std::sqrt(pow((0.03/HitDetect.Length()), 2.0)+pow((timeRes/tof), 2.0));
						Eeject += rndgauss0(energyRes*Eeject); // Smear energy due to VANDLE resolution
					}
				
					// Main output string
					// X(m) Y(m) Z(m) LabTheta(deg) LabPhi(deg) EjectileE(MeV) EjectileToF(ns) Bar# Face# HitX(m) HitY(m) HitZ(m)
					if(Eeject >= 0.1 && Eeject <= 3.0){
						VIKARout << HitDetect.axis[0] << "\t" << HitDetect.axis[1] << "\t" << HitDetect.axis[2] << "\t"; // Hit X, Y, and Z (in meters)
						VIKARout << HitDetectSphere.axis[1]*rad2deg << "\t" << HitDetectSphere.axis[2]*rad2deg << "\t" << Eeject << "\t" << tof*(1E9) << "\t"; // Energy, ToF, and ejectile angles
						VIKARout << bar << "\t" << face << "\t" << hit_x << "\t" << hit_y << "\t" << hit_z << std::endl; // Individual bar data
						if(WriteRecoil){ RecoilOut << RecoilSphere.axis[1] << "\t" << RecoilSphere.axis[2] << "\t" << Erecoil << "\n"; }
			
						Ndetected++;
						if(!flag){ flag = true; }
					}
				}

				break;
			}
		}

	} // Main simulation loop
	// ==  ==  ==  ==  ==  ==  == 
 
	// Information output and cleanup
	std::cout << "\n ------------- Simulation Complete --------------\n";
	std::cout << " Simulation Time: " << (float)(clock()-timer)/CLOCKS_PER_SEC << " seconds\n"; 
	if(!PerfectDet){
		std::cout << " Geometric Efficiency: " << NbarHit*100.0/Nreactions << "%\n";
		std::cout << " Detection Efficiency: " << Ndetected*100.0/Nreactions << "%\n"; 
	}
	else{ std::cout << " Detection Efficiency: " << Ndetected*100.0/Nreactions << "%\n"; }
	if(SupplyRates){ std::cout << " Beam Time: " << Nsimulated/BeamRate << " seconds\n"; }
	
	std::cout << "  Wrote file " << output_fname_prefix+"_main.dat\n"; VIKARout.close();
	if(WriteRecoil){ std::cout << "  Wrote file " << output_fname_prefix+"_recoil.dat\n"; RecoilOut.close(); }
	if(WriteDebug){ std::cout << "  Wrote file " << output_fname_prefix+"_debug.dat\n"; DebugOut.close(); }
} 
