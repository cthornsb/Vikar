// vikar.h
// Converted by FortranConvert v0.1
// Wed Feb 12 19:55:17 2014

#include <fstream>
#include <iostream>
#include <time.h>

#include "vikar_core.h"

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
		
	std::ofstream file22("Info.dat"); // (CORY) For debugging
				
	// Arrays for spectra
	int testArray[360][360]; 
	
	// temporary variables
	double thetaEjectdE, PhiEjectdE; 
	double thetaRecoildE, PhiRecoildE; 
	
	// Energy loss varibales
	bool SRIM_conv; 
	unsigned short NtargElements, beamTables, ejectTables_targ, recoilTables_targ; 
	unsigned short N_SRIMpoints, ejectTables_det, recoilTables_det, **targComp; 
	double tgtdens, energy, dedxmg; 
	double TgtIonPot, rangemg; 
	double thickness, Zdepth, targ_depth, targ_angle; 
	double range_beam;
	double det_E_eject, det_dE_eject; 
	double det_E_recoil, det_dE_recoil; 
	double Ereact, Eemerge_eject, Eemerge_recoil; 
	double AveTargZ, AveTargA, totalTarg; 
	double targRadL, DetRadL, dummy; 
	double density_det, conv_det, Adet, Zdet; // details of detector material for loss/straggling calcs
	
	// Energy loss arrays
	double *beam_Erange, *beam_E, *eject_Erange, *eject_E, *eject_Erange_det, *eject_E_det;
	double *recoil_Erange, *recoil_E, *recoil_Erange_det, *recoil_E_det, *beam_long, *beam_lat;
	double *eject_long, *eject_lat, *eject_long_det, *eject_lat_det, *recoil_long, *recoil_lat;
	double *recoil_long_det, *recoil_lat_det, *SRIM_dedx;
	
	// Physics Variables
	unsigned short body, NRecoilStates, RecoilState, ADists; // Number of body reaction, no. states in recoil populated
	std::string *AngDist_fname; 
	unsigned short *ADist_type, *NAngPoints;
	double *ExRecoilStates, gsQvalue; 
	double **DistAng, **DistInt, *maxInt; 
	double *DistAng_temp, *DistInt_temp, *totXsect; 
	double Ebeam, Ebeam0, Erecoil, Eeject, Abeam, Zbeam, Atarg;
	double Arecoil, Zrecoil, Aeject, Zeject, Q, ran_Ex_recoil;
	double Ztarg;
	
	// Beam variables
	double thetaBeam, phiBeam; // Incident beam particle direction
	double thetaReact, phiReact; // Direction of beam particle at reaction (after straggling)
	double beamspot, beamEspread; // Beamspot size FWHM (mm)
	
	// Reaction product variables
	double thetaRecoil, phiRecoil; // Direction of recoils after reaction
	double thetaEject, phiEject; // Direction of ejectiles after reaction
	double ExRecoil, ExEject; // Excitation of the reaction products
	
	//double lost_E_eject; // Energy deposited by the ejectile in the E detector
	double ThetaEjectE, PhiEjectE; // Apparent direction of the ejectile after dE straggling
	double thetaRecoilE, PhiRecoilE; // Apparent direction of the recoil after dE straggling
	double dE_final_Theta_eject, E_final_Theta_eject; // Final (smeared) theta for the ejectile from the dE and E
	double dE_final_Phi_eject, E_final_Phi_eject; // Final (smeared) phi for the ejectile from the dE and E
	double dE_final_eject, E_final_eject; // Final (smeared) energies for the ejectile from the dE and E
	
	double dE_final_Theta_recoil, E_final_Theta_recoil; // Final (smeared) theta for the recoil from the dE and E
	double dE_final_Phi_recoil, E_final_Phi_recoil; // Final (smeared) phi for the recoil from the dE and E
	double dE_final_recoil, E_final_recoil; // Final (smeared) energies for the recoil from the dE and E
	
	// Detector variables
	std::string DetSetCyl_fName, DetSetPlan_fName, DetSetAnn_fName; 
	bool Det_eject, Det_recoil, coincidence;
	bool *xPlane, *yPlane; 
	
	unsigned short detNo, stripNo, DetSetup; 
	unsigned short perfectDet; 
	unsigned short EjectHit, Eject_hitGeo, RecoilHit, Recoil_hitGeo; 
	unsigned short StripHitCyl_eject, StripHitPlan_eject, StripHitAnn_eject; 
	unsigned short StripHitCyl_recoil, StripHitPlan_recoil, StripHitAnn_recoil; 
	unsigned short DetHitCyl_eject, DetHitPlan_eject, DetHitAnn_eject; // Need these?
	unsigned short DetHitCyl_recoil, DetHitPlan_recoil, DetHitAnn_recoil; // Need these?
	unsigned short NDet_cyl, Ndet_plan, Ndet_ann; 
	unsigned short *Nstrips_cyl, *Nstrips_plan, *Nstrips_ann;
	double *detLength_cyl, *detWidth_cyl, *DetRad_cyl, *DetZ_cyl, *DetPhi_cyl; 
	double **DetThetaMax, **DetThetaMin, **DetPhiMax, **DetPhiMin; 
	double *dEthick_cyl, *Ethick_cyl, *Eres_dE_cyl, *Pres_dE_cyl, *Eres_E_cyl; 
	double *Pres_E_cyl, *Pres_dE_en_cyl, *Pres_E_en_cyl; // energy at which Pres is measured
	
	double *DetZ_plan, *DetX_plan, *DetY_plan, *detLength_plan, *detWidth_plan; 
	double *dEthick_plan, *Ethick_plan, *Eres_dE_plan, *Pres_dE_plan; 
	double *Pres_dE_en_plan, *Pres_E_en_plan, *Eres_E_plan, *Pres_E_plan; // energy at which Pres is measured
	double **DetXMax, **DetXMin, **DetYMax, **DetYMin, *XMax, *XMin, *YMax, *YMin; 
	
	double *DetZ_ann, *DetInner_ann, *DetOuter_ann, *detPhi_ann; //*detdPhi_ann; 
	double *dEthick_ann, *Ethick_ann, *Eres_dE_ann, *Eres_E_ann; 
	double **DetThetaMax_ann, **DetThetaMin_ann, **DetPhiMax_ann, **DetPhiMin_ann;
	
	// Input/output variabless
	unsigned short num_beam_E, num_eject_E, num_recoil_E;
	unsigned short idummy, num_eject_det_E, num_recoil_det_E;
	int Ndetected, Nwanted, Nsimulated;
	bool write_eject, write_recoil; 
	bool write_all_eject, write_all_recoil; 
	unsigned short n, p, i, j;
	clock_t timer;

	// Error reporting variables
	unsigned short version, version_check; 
	float version_read; 
	bool SRIM_stat, stuck; 
	std::string SRIM_fName, SRIM_fName_beam, SRIM_fName_targ_eject, SRIM_fName_targ_recoil;
	std::string SRIM_fName_det_eject, SRIM_fName_det_recoil;
	//char SRIM_fName[30]; 

	bool wfile = false;
	
	// Zero the test output array
	for(i = 0; i <= 360; i++){ 
		for(j = 0; j <= 360; j++){ 
			testArray[i][j] = 0; 
		} 
	} 

	n = 1; 
	p = 1; 

	//------------------------------------------------------------------------
	//
	// Start of user input code
	//
	//------------------------------------------------------------------------

	version = 312; 

	std::cout << "\n####       #### ######## ####     ###         ###       #########\n"; 
	std::cout << " ##         ##     ##     ##     ##          ## ##       ##     ##\n";
	std::cout << "  ##       ##      ##     ##    ##          ##   ##      ##     ##\n";
	std::cout << "  ##       ##      ##     ##   ##           ##   ##      ##   ##\n";
	std::cout << "   ##     ##       ##     #####            #########     #####\n";
	std::cout << "   ##     ##       ##     ##   ##          ##     ##     ##   ##\n";
	std::cout << "    ##   ##        ##     ##    ##        ##       ##    ##    ##\n";
	std::cout << "    ##   ##        ##     ##     ##       ##       ##    ##     ##\n";
	std::cout << "     ## ##         ##     ##      ##     ##         ##   ##      ##\n";
	std::cout << "      ###       ######## ####      ###  ####       #### ####      ###\n";

	std::cout << "\n VIKAR 3.12\n"; 
	std::cout << " ==  ==  ==  ==  == \n\n"; 
	
	std::cout << " Welcome to VIKAR << the Virtual Instrumentation\n"; 
	std::cout << " for Kinematics And Reactions program\n\n"; 

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
		
		unsigned short count = 0;
		std::string input;
		std::cout << "\n Reading from file " << argv[1] << std::endl;
		while(!input_file.eof()){
			input_file >> input;
			if(count == 0){  }
			else if(count == 1){ 
				version = short(atof(input.c_str())*100); 
				std::cout << "  Version: " << input << std::endl;
			}
			else if(count == 2){ 
				body = atos(input.c_str()); 
				std::cout << "  No. Body: " << body << std::endl;
			}
			else if(count == 3){ 
				Zbeam = atof(input.c_str());  
				std::cout << "  Beam-Z: " << Zbeam << std::endl;
			}
			else if(count == 4){ 
				Abeam = atof(input.c_str());  
				std::cout << "  Beam-A: " << Abeam << std::endl;
			}
			else if(count == 5){ 
				Ztarg = atof(input.c_str());  
				std::cout << "  Target-Z: " << Ztarg << std::endl;
			}
			else if(count == 6){ 
				Atarg = atof(input.c_str());  
				std::cout << "  Target-A: " << Atarg << std::endl;
			}
			else if(count == 7){ 
				Zeject = atof(input.c_str());
				Zrecoil = Zbeam + Ztarg - Zeject;  
				std::cout << "  Ejectile-Z: " << Zeject << std::endl;
			}
			else if(count == 8){ 
				Aeject = atof(input.c_str()); 
				Arecoil = Abeam + Atarg - Aeject; 
				std::cout << "  Ejectile-A: " << Aeject << std::endl;
				std::cout << "  Recoil-Z: " << Zrecoil << std::endl;
				std::cout << "  Recoil-A: " << Arecoil << std::endl;
			}
			else if(count == 9){ 
				Ebeam = atof(input.c_str());  
				std::cout << "  Beam Energy: " << Ebeam << " MeV\n";
			}
			else if(count == 10){ 
				beamspot = atof(input.c_str());  
				std::cout << "  Beam Spot Size: " << beamspot << " mm\n";
				beamspot *= 0.1;
			}
			else if(count == 11){ 
				beamEspread = atof(input.c_str());  
				std::cout << "  Beam Spread: " << beamEspread << " MeV\n";
			}
			else if(count == 12){ 
				gsQvalue = atof(input.c_str());  
				std::cout << "  G.S. Q-Value: " << gsQvalue << " MeV\n";
			}
			else if(count == 13){ 
				NRecoilStates = atos(input.c_str());
				std::cout << "  No. Recoil States: " << NRecoilStates << std::endl;
				ExRecoilStates = new double[NRecoilStates];
				totXsect = new double[NRecoilStates];
				for(i = 0; i < NRecoilStates; i++){
					input_file >> input;
					ExRecoilStates[i] = atof(input.c_str());
					std::cout << "   Recoil State " << i+1 << ": " << ExRecoilStates[i] << " MeV\n";
					totXsect[i] = 0.0;
				}
			}
			else if(count == 14){ 
				ADists = atos(input.c_str());
				std::cout << "  Supply Angular Distributions: ";
				ADist_type = new unsigned short[NRecoilStates];
				AngDist_fname = new std::string[NRecoilStates];
				NAngPoints = new unsigned short[NRecoilStates];
				maxInt = new double[NRecoilStates];
				if(ADists == 1){
					std::cout << "Yes\n";
					for(i = 0; i < NRecoilStates; i++){
						input_file >> input;
						AngDist_fname[i] = input;
						std::cout << "   Distribution for state " << i+1 << ": " << AngDist_fname[i] << std::endl;
						AngDist_read(AngDist_fname[i].c_str(),NAngPoints[i],DistAng_temp,DistInt_temp,maxInt[i]); 
						totXsect[i] = totXsect[i-1]+maxInt[i]; 
						if(NAngPoints[i] > 0){
							ADist_type[i] = 1; // set to use this distribution
							for(j = 0; j < NAngPoints[i]; j++){
								DistAng[i][j] = DistAng_temp[j]; 
								DistInt[i][j] = DistInt_temp[j]; 
							} 
						} 
						else{  ADist_type[i] = 0; } // force to isotropic as distribution non-existant!
					}
				}
				else{
					std::cout << "No\n";
					for(i = 0; i < NRecoilStates; i++){
						ADist_type[i] = 0;
						AngDist_fname[i] = "";
						NAngPoints[i] = 0;
						maxInt[i] = 0;
					}
				}
			}
			else if(count == 15){ 
				NtargElements = atos(input.c_str());
				std::cout << "  No. Target Elements: " << NtargElements << std::endl;
				targComp = new unsigned short *[NtargElements];
				for(i = 0; i < NtargElements; i++){
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
				for(i = 0; i < NtargElements; i++){
					AveTargZ = AveTargZ + targComp[i][0]*targComp[i][2]; 
					AveTargA = AveTargA + targComp[i][1]*targComp[i][2]; 
					totalTarg = totalTarg + targComp[i][2]; 
					dummy = radlength(targComp[i][1],targComp[i][0]); 
					targRadL = targRadL + ((targComp[i][1]*targComp[i][2]/AveTargA)/dummy); 
				} 
	
				targRadL = 1.0/targRadL; 
				AveTargZ = AveTargZ/totalTarg; 
				AveTargA = AveTargA/totalTarg;
			}
			else if(count == 16){ 
				thickness = atof(input.c_str());  
				std::cout << "  Target Thickness: " << thickness << " mg/cm^2\n";
			}
			else if(count == 17){ 
				tgtdens = atof(input.c_str());  
				std::cout << "  Target Density: " << tgtdens << " g/cm^3\n";
			}
			else if(count == 18){ 
				targ_angle = atof(input.c_str());  
				std::cout << "  Target Angle: " << targ_angle << " degrees\n";
				targ_angle *= deg2rad;
			}
			else if(count == 19){ 
				beamTables = atos(input.c_str()); 
				std::cout << "  Beam Range Table: ";
				if(beamTables == 1){ 
					input_file >> input;
					SRIM_fName_beam = input; 
					std::cout << SRIM_fName_beam << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 20){ 
				ejectTables_targ = atos(input.c_str());
				std::cout << "  Ejectile Range Table (Target): ";
				if(ejectTables_targ == 1){
					input_file >> input;
					SRIM_fName_targ_eject = input;
					std::cout << SRIM_fName_targ_eject << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 21){ 
				recoilTables_targ = atos(input.c_str());
				std::cout << "  Residual Range Table (Target): ";
				if(recoilTables_targ == 1){
					input_file >> input;
					SRIM_fName_targ_recoil = input;
					std::cout << SRIM_fName_targ_recoil << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 22){ 
				ejectTables_det = atos(input.c_str());
				std::cout << "  Ejectile Range Table (Detector): ";
				if(ejectTables_det == 1){
					input_file >> input;
					SRIM_fName_det_eject = input;
					std::cout << SRIM_fName_det_eject << std::endl;
				}
				else{ std::cout << "No\n"; }
			}	
			else if(count == 23){ 
				recoilTables_det = atos(input.c_str());
				std::cout << "  Residual Range Table (Detector): ";
				if(recoilTables_det == 1){
					input_file >> input;
					SRIM_fName_det_recoil = input;
					std::cout << SRIM_fName_det_recoil << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 24){ 
				perfectDet = atos(input.c_str()); 
				std::cout << "  Perfect Detector: ";
				if(perfectDet == 1){ std::cout << "Yes\n"; }
				else{ std::cout << "No\n"; }
			}
			else if(count == 25){ 
				DetSetup = atos(input.c_str());
				std::cout << "  Detector Setup File: ";
				if(DetSetup == 1){
					std::cout << "Yes\n";
					input_file >> input; std::cout << "   Cylindrical: " << input << std::endl;
					input_file >> DetSetPlan_fName; std::cout << "   Planar: " << DetSetPlan_fName << std::endl;
					input_file >> input; std::cout << "   Annular: " << input << std::endl;
				}
				else{ std::cout << "No\n"; }
			}
			else if(count == 26){ 
				Nwanted = atol(input.c_str());
				std::cout << "  Desired Detections: " << Nwanted << std::endl; 
			}
			else if(count == 27){ 
				idummy = atos(input.c_str());
				std::cout << "  Detect Ejectiles: "; 
				if(idummy == 1){ std::cout << "Yes\n"; Det_eject = true; }
				else{ std::cout << "No\n"; Det_eject = false; }
			}
			else if(count == 28){ 
				idummy = atos(input.c_str()); 
				std::cout << "  Detect Recoils: "; 
				if(idummy == 1){ std::cout << "Yes\n"; Det_recoil = true; }
				else{ std::cout << "No\n"; Det_recoil = false; }
				
				if(Det_eject && Det_recoil){ 
					input_file >> input; 
					idummy = atos(input.c_str());
					std::cout << "  Require Coincidence: ";
					if(idummy == 1){ std::cout << "Yes\n"; coincidence = true; }
					else{ std::cout << "No\n"; coincidence = false; }
				}
			}
			else if(count == 29){ 
				idummy = atos(input.c_str());
				std::cout << "  Write All Ejectiles: "; 
				if(idummy == 1){ std::cout << "Yes\n"; write_all_eject = true; }
				else{ std::cout << "No\n"; write_all_eject = false; }					
			}
			else if(count == 30){ 
				idummy = atos(input.c_str()); 
				std::cout << "  Write All Residuals: "; 
				if(idummy == 1){ std::cout << "Yes\n"; write_all_recoil = true; }
				else{ std::cout << "No\n"; write_all_recoil = false; }
			}
			
			count++;
		}
		
		input_file.close();
		if(count < 30){
			std::cout << " Error: The input file is invalid\n";
			return 1;
		}
	}
	else{
		std::cout << " Error: Missing required variable\n";
		return 1;
	}
	
	std::cout << "\n ==  ==  ==  ==  == \n\n";
	
	if(beamTables == 1){ 		
		std::cout << " Reading in range file for beam particles...\n"; 
		std::cout << N_SRIMpoints << std::endl;
		SRIMread(SRIM_fName_beam,SRIM_stat,N_SRIMpoints,beam_E,SRIM_dedx,beam_Erange,beam_long,beam_lat,true); 
	
		// Check whether SRIM_stat returned a good error code
		if(!SRIM_stat){
			std::cout << " Error: SRIMread() returned an error code for beamTables!\n";
			return 1;
		} 
		
		// Determine the ranges that straddle the beam energy
		// and calculate the range at the specific beam energy
		i = 0; 
		while((Ebeam+10*beamEspread)-beam_E[i] > 0){ i = i+1; } 
		range_beam = linear(beam_E[i-1],beam_E[i],beam_Erange[i-1],beam_Erange[i],Ebeam); 
	}
	else{ 
		// If internal range tables are to be used for beam particles
		// Calculate the stopping power table for the beam particles in the target
		std::cout << " Calling internal routine for beam particle ranges...\n";
		num_beam_E = short(Ebeam*1.0) + short(10*beamEspread) + 10;
		beam_Erange = new double[num_beam_E];
		beam_E = new double[num_beam_E];
		for(i = 0; i < num_beam_E; i++){
			energy = (i+1)/1.0; 
          		ncdedx(tgtdens,AveTargA,AveTargZ,Abeam,Zbeam,energy,dedxmg,TgtIonPot,rangemg);
          		beam_Erange[i] = rangemg;
          		beam_E[i] = energy;
		} 
	}
	
	if(ejectTables_targ == 1){
		std::cout << " Reading in range file for ejectiles...\n"; 
		SRIMread(SRIM_fName,SRIM_stat,N_SRIMpoints,eject_E,SRIM_dedx,eject_Erange,eject_long,eject_lat,SRIM_conv); 
	
		// Check whether SRIM_stat returned a good error code
		if(!SRIM_stat){
			std::cout << " Error: SRIMread() returned an error code for ejectTables_targ!\n";
			return 1;
		} 
	}
	else{ 
		// Calculate the stopping power table for the ejectiles in the target
		std::cout << " Calling internal routine for ejectile ranges...\n";
		num_eject_E = short((Ebeam + 10*beamEspread + 10)*10.0);
		eject_Erange = new double[num_eject_E];
		eject_E = new double[num_eject_E];
		for(i = 0; i < num_eject_E; i++){
          		energy = (i+1)/10.0;
			ncdedx(tgtdens,AveTargA,AveTargZ,Aeject,Zeject,energy,dedxmg,TgtIonPot,rangemg);
			eject_Erange[i] = rangemg; 
			eject_E[i] = energy; 
		} 
	}

	if(recoilTables_targ == 1){
		std::cout << " Reading in range file for recoils...\n"; 
		SRIMread(SRIM_fName,SRIM_stat,N_SRIMpoints,recoil_E,SRIM_dedx,recoil_Erange,recoil_long,recoil_lat,SRIM_conv); 
	
		// Check whether SRIM_stat returned a good error code
		if(!SRIM_stat){
			std::cout << " Error: SRIMread() returned an error code for recoilTables_targ!\n";
			return 1;
		} 
	}
	else{ 
		// Calculate the stopping power table for the recoils in the target
		std::cout << " Calling internal routine for recoil ranges...\n";
		num_recoil_E = short(Ebeam*1.0) + short(10*beamEspread) + 100;
		recoil_Erange = new double[num_recoil_E];
		recoil_E = new double[num_recoil_E];
		for(i = 0; i < num_recoil_E; i++){ 
			energy = (i+1)/1.0; 
			ncdedx(tgtdens,AveTargA,AveTargZ,Arecoil,Zrecoil,energy,dedxmg,TgtIonPot,rangemg); 
			recoil_Erange[i] = rangemg; 
			recoil_E[i] = energy; 
		} 
	}

	// Set detector material to Si. Add user selection at some point
	Adet = 28; 
	Zdet = 14; 
	density_det = 2.3212; 
	density_det = density_det*1000.0; // Convert density to mg/cm^3
	conv_det = 1.0e-4*density_det; 
	DetRadL = radlength(Adet,Zdet); 
	
	// If external range tables are to be used for ejectiles
	if(ejectTables_det == 1){
		std::cout << " Reading in range file for ejectiles..."; 
		SRIMread(SRIM_fName,SRIM_stat,N_SRIMpoints,eject_E_det,SRIM_dedx,eject_Erange_det,eject_long_det,eject_lat_det,SRIM_conv); 
	
		// Check whether SRIM_stat returned a good error code
		if(!SRIM_stat){
			std::cout << " Error: SRIMread() returned an error code for ejectTables_det!\n";
			return 1;
		} 
	}
	else{ 
		// Calculate the stopping power table for the ejectiles in the target
		std::cout << " Calling internal routine for ejectile ranges...\n"; 
		num_eject_det_E = short((Ebeam + 10*beamEspread + 10)*10.0);
		eject_Erange_det = new double[num_eject_det_E];
		eject_E_det = new double[num_eject_det_E];
		for(i = 0; i < num_eject_det_E; i++){ 
			energy = (i+1)/10.0; 
			std::cout << tgtdens << "," << Adet << "," << Zdet << "," << Aeject << "," << Zeject << "," << energy << "," << dedxmg << std::endl;
          		ncdedx(tgtdens,Adet,Zdet,Aeject,Zeject,energy,dedxmg,TgtIonPot,rangemg);
          		eject_Erange_det[i] = rangemg/conv_det; 
			eject_E_det[i] = energy; 
		} 
	}

	// If external range tables are to be used for recoils
	if(recoilTables_det == 1){
		std::cout << " Reading in range file for recoils...\n"; 
		SRIMread(SRIM_fName,SRIM_stat,N_SRIMpoints,recoil_E_det,SRIM_dedx,recoil_Erange_det,recoil_long_det,recoil_lat_det,SRIM_conv);
	
		// Check whether SRIM_stat returned a good error code
		if(!SRIM_stat){
			std::cout << " Error: SRIMread() returned an error code for recoilTables_det!\n";
			return 1;
		} 
	}
	else{ 
		// Calculate the stopping power table for the recoils in the target
		std::cout << " Calling internal routine for recoil ranges...\n"; 
		num_recoil_det_E = short((Ebeam + 10*beamEspread + 10)*10.0);
		recoil_Erange_det = new double[num_recoil_det_E];
		recoil_E_det = new double[num_recoil_det_E];
		for(i = 0; i < num_recoil_det_E; i++){
			energy = (i+1)/10.0; 
         		ncdedx(tgtdens,Adet,Zdet,Arecoil,Zrecoil,energy,dedxmg,TgtIonPot,rangemg);
			recoil_Erange_det[i] = rangemg/conv_det; 
			recoil_E_det[i] = energy; 
		} 
	}
	
	std::cout << "\n ==  ==  ==  ==  == \n\n";

	if(perfectDet == 1 && DetSetup == 1){
		/*std::cout << " Reading in cylindrical detector setup file...\n"; 
		DetSet_cyl_read(DetSetCyl_fName,NDet_cyl,Nstrips_cyl,detLength_cyl,DetWidth_cyl,DetRad_cyl,DetZ_cyl,DetPhi_cyl,dEthick_cyl,Ethick_cyl,Eres_dE_cyl,Pres_dE_cyl,Eres_E_cyl,Pres_E_cyl,Pres_dE_en_cyl,Pres_E_en_cyl); 

		// Report on how many cylindrical detectors were read in
		std::cout << "  Found " << NDet_cyl << " cylindrical detectors in file " << DetSetCyl_fName << std::endl;*/

		Ndet_plan = GetLines(DetSetPlan_fName.c_str());
		detLength_plan = new double[Ndet_plan];
		detWidth_plan = new double[Ndet_plan];
		Nstrips_plan = new unsigned short[Ndet_plan];
		DetZ_plan = new double[Ndet_plan];
		DetX_plan = new double[Ndet_plan];
		DetY_plan = new double[Ndet_plan];
		dEthick_plan = new double[Ndet_plan];
		Ethick_plan = new double[Ndet_plan];
		Eres_dE_plan = new double[Ndet_plan];
		Pres_dE_plan = new double[Ndet_plan];
		Eres_E_plan = new double[Ndet_plan];
		Pres_E_plan = new double[Ndet_plan];
		Pres_dE_en_plan = new double[Ndet_plan];
		Pres_E_en_plan = new double[Ndet_plan];
		xPlane = new bool[Ndet_plan];
		yPlane = new bool[Ndet_plan];
		
		std::cout << " Reading in planar detector setup file...\n";
		DetSet_plan_read(DetSetPlan_fName.c_str(),Ndet_plan,Nstrips_plan,detLength_plan,detWidth_plan,DetZ_plan,DetX_plan,DetY_plan,xPlane,yPlane,
				 dEthick_plan,Ethick_plan,Eres_dE_plan,Pres_dE_plan,Eres_E_plan,Pres_E_plan,Pres_dE_en_plan,Pres_E_en_plan); 

		// Report on how many detectors were read in
		std::cout << " Found " << Ndet_plan << " planar detectors in file " << DetSetPlan_fName << std::endl;

		/*std::cout << " Reading in annular detector setup file...\n"; 
		DetSet_ann_read(DetSetAnn_fName,Ndet_ann,Nstrips_ann,detInner_ann,DetOuter_ann,DetdPhi_ann,DetPhi_ann,DetZ_ann,dEthick_ann,Ethick_ann,Eres_dE_ann,Eres_E_ann); 

		// Report on how many detectors were read in
		std::cout << " Found " << Ndet_ann << " annular detectors in file " << DetSetAnn_fName << std::endl;*/

		// Check there's at least 1 detector!
		if((NDet_cyl + Ndet_plan + Ndet_ann) < 1){
			std::cout << " Error: Found no detectors. Check that the filenames are correct\n"; 
			return 1;
		}
	}

	// Detector Setup format lines
	// XXXXXX Put this in a subroutine, then can it as needed????
	// Determine the angular range covered by each strip of each cylindrical detector
	/*for(detNo = 0; detNo < NDet_cyl; detNo++){
		DetAng_cyl(detLength_cyl[detNo],detWidth_cyl[detNo],Nstrips_cyl[detNo],DetRad_cyl[detNo],DetZ_cyl[detNo],DetPhi_cyl[detNo],thetaMax,thetaMin,PhiMax,PhiMin); 
		for(stripNo = 0; stripNo < Nstrips_cyl[detNo]; stripNo++){
			DetThetaMax[detNo][stripNo] = thetaMax[stripNo]; 
			DetThetaMin[detNo][stripNo] = thetaMin[stripNo]; 
			DetPhiMax[detNo][stripNo] = PhiMax[stripNo]; 
			DetPhiMin[detNo][stripNo] = PhiMin[stripNo]; 
			// "'---------------'";
			// "'Detector ',detNo,' Strip ',stripNo";
			// "'DetThetaMax = ',DetThetaMax[detNo][stripNo]*rad2deg";
			// "'DetThetaMin = ', DetThetaMin[detNo][stripNo]*rad2deg";
			// "'DetPhiMax = ',DetPhiMax[detNo][stripNo]*rad2deg";
			// "'DetPhiMin = ',DetPhiMin[detNo][stripNo]*rad2deg";
		} 
	} // Over detNo*/
	
	DetXMax = new double *[Ndet_plan];
	DetXMin = new double *[Ndet_plan];
	DetYMax = new double *[Ndet_plan];
	DetYMin = new double *[Ndet_plan];
	for(detNo = 0; detNo < Ndet_plan; detNo++){
		DetXMax[detNo] = new double[Nstrips_plan[detNo]];
		DetXMin[detNo] = new double[Nstrips_plan[detNo]];
		DetYMax[detNo] = new double[Nstrips_plan[detNo]];
		DetYMin[detNo] = new double[Nstrips_plan[detNo]];
		DetAng_plan(detLength_plan[detNo],detWidth_plan[detNo],Nstrips_plan[detNo],DetZ_plan[detNo],DetX_plan[detNo],DetY_plan[detNo],
			    xPlane[detNo],yPlane[detNo],DetXMax[detNo],DetXMin[detNo],DetYMax[detNo],DetYMin[detNo]); 
	} 

	/*for(detNo = 0; detNo < ndet_ann; detNo++){
		DetAng_ann(detInner_ann[detNo],detOuter_ann[detNo],detdPhi_ann[detNo],detPhi_ann[detNo],Nstrips_ann[detNo],DetZ_ann[detNo],thetaMax,thetaMin,PhiMax,PhiMin); 
		for(stripNo = 0; stripNo < Nstrips_ann[detNo]; stripNo++){
			DetThetaMax_ann[detNo][stripNo] = thetaMax[stripNo]; 
			DetThetaMin_ann[detNo][stripNo] = thetaMin[stripNo]; 
			DetPhiMax_ann[detNo][stripNo] = PhiMax[stripNo]; 
			DetPhiMin_ann[detNo][stripNo] = PhiMin[stripNo]; 
			// "'---------------'";
			// "'Detector ',detNo,' Strip ',stripNo";
			// "'DetThetaMax = ',DetThetaMax_ann[detNo][stripNo]*rad2deg";
			// "'DetThetaMin = ',DetThetaMin_ann[detNo][stripNo]*rad2deg";
			// "'DetPhiMax = ',DetPhiMax_ann[detNo][stripNo]*rad2deg";
			// "'DetPhiMin = ',DetPhiMin_ann[detNo][stripNo]*rad2deg";
		} 
	} // Over detNo*/

	std::cout << "\n ==  ==  ==  ==  == \n\n";

	//---------------------------------------------------------------------------
	// End of Input Section
	//---------------------------------------------------------------------------
	//e = dtime(t); // incompatible

	Ndetected = 0; 
	Nsimulated = 0;

	// Open some output files
	std::ofstream file11("unsmeared_all.dat"); 
	std::ofstream file12("unsmeared.dat"); 
	std::ofstream file13("angles.dat"); 
	std::ofstream file14("Etot_thetaE.dat"); 
	std::ofstream file15("Etot_thetadE.dat"); 
	std::ofstream file16("dE_thetaE.dat"); 
	std::ofstream file17("dE_thetadE.dat"); 
	std::ofstream file18("E_thetaE.dat"); 
	std::ofstream file19("E_thetadE.dat"); 

	std::ofstream file20("PI-f.dat"); 
	std::ofstream file21("PI-b.dat"); 
	//      open(unit = 22,file = 'out_recoil_PI-f.dat',status = 'unknown')
	//      open(unit = 23,file = 'out_recoil_PI-b.dat',status = 'unknown')

	//      open(unit = 41,file = 'out_recoil.dat',status = 'unknown')
	//      open(unit = 41,file = 'out_recoil_unsmeared_all.dat',status = 'unknown')
	//      open(unit = 42,file = 'out_recoil_unsmeared.dat',status = 'unknown')
	//      open(unit = 43,file = 'out_recoil_angles.dat',status = 'unknown')
	//      open(unit = 44,file = 'out_recoil_Etot_thetaE.dat',status = 'unknown')
	//      open(unit = 45,file = 'out_recoil_Etot_thetadE.dat',status = 'unknown')
	//      open(unit = 46,file = 'out_recoil_dE_thetaE.dat',status = 'unknown')
	//      open(unit = 47,file = 'out_recoil_dE_thetadE.dat',status = 'unknown')
	//      open(unit = 48,file = 'out_recoil_E_thetaE.dat',status = 'unknown')
	//      open(unit = 49,file = 'out_recoil_E_thetadE.dat',status = 'unknown')

	//      open(unit = 57,file = 'debug.dat',status = 'unknown')
	//      open(unit = 58,file = 'debug2.dat',status = 'unknown')

	// Begin the simulation
	std::cout << " ---------- Simulation Setup Complete -----------\n"; 
	std::cout << "\n Beginning simulating a bunch of events....\n"; 

	Ebeam0 = Ebeam; 

	//---------------------------------------------------------------------------
	// The Event Loop
	// ==  ==  ==  ==  ==  ==  == 
	// (Just to make it obvious)
	//---------------------------------------------------------------------------

	int chunk = Nwanted/10;
	float totTime = 0.0;
	char counter = 1;
	bool flag = false;
	timer = clock();

	while(Ndetected < Nwanted){
		// ****************Time Estimate**************
		if(flag && (Ndetected % chunk == 0)){
			flag = false;
			totTime = (float)(clock()-timer)/CLOCKS_PER_SEC;
			std::cout << "\n ------------------------------------------------\n"; 
			std::cout << " Number of particles Simulated: " << Nsimulated << std::endl; 
			std::cout << " Number of particles Detected: " << Ndetected << std::endl; 
			std::cout << " " << Ndetected*100.0/Nwanted << "% of simulation complete...\n"; 
			std::cout << "  Detection Efficiency: " << Ndetected*100.0/Nsimulated << "%\n"; 
			std::cout << "  Time taken = " << totTime << " seconds\n";
			std::cout << "  Time reamining: " << (totTime/counter)*(10-counter) << " seconds\n";
			counter++; 
		}

		// Recoil Hit????

		// Set the beam direction to (0,0,1) i.e. along the z axis

		// Set the incident direction for the beam
		//	thetaBeam = frand()/180.0*pi*10.0
		thetaBeam = rndgauss0(0.1*deg2rad); 
		phiBeam = frand()*2.0*pi; 

		// Calculate the depth in the target for the reaction
		targ_depth = thickness *frand(); 

		// Calculate thickness traversed by the beam particle (Zdepth),
		// dependent on the target angle and beam particle direction.
		// (use (thickness-targ_depth) as subroutine written for recoil/ejectile)
		dummy = thickness-targ_depth; 
		targ_thick(thetaBeam,phiBeam,thickness,dummy,targ_angle,Zdepth); 

		// Calculate the beam particle energy, varied with energy spread
		Ebeam = Ebeam0 + rndgauss0(beamEspread); 

		// Calculate the beam particle range in the target
		//range_beam = linear(beam_E[short(Ebeam)],beam_E[short(Ebeam+1.0)],beam_Erange[short(Ebeam)],beam_Erange[short(Ebeam+1.0)],Ebeam); 
		range_beam = linear(beam_E[short(Ebeam)-1],beam_E[short(Ebeam)],beam_Erange[short(Ebeam)-1],beam_Erange[short(Ebeam)],Ebeam); 

		// Code for calculating energy loss of particle...
		// Find the ranges for energies straddling the beam energy
		//i = 1; 
		//while(beam_Erange[i] < (range_beam-Zdepth)){ i = i + 1; }
		for(i = 1; i < num_beam_E; i++){
			if(beam_Erange[i] < (range_beam-Zdepth)){ break; }
		}

		// Calculate the new energy
		Ereact = linear(beam_Erange[i-1],beam_Erange[i],beam_E[i-1],beam_E[i],(range_beam-Zdepth)); 
		//std::cout << i << " " << beam_Erange[i-1] << " " << beam_Erange[i] << " " << beam_E[i-1] << " " << beam_E[i] << " " << range_beam << " " << Zdepth << std::endl; // CORY!

		// Ereact = Ebeam // Uncomment to remove beam particle energy loss effects
		//Eloss = Ebeam - Ereact; 

		// Determine the angle of the beam particle's trajectory at the
		// interaction point, due to angular straggling and the incident
		// trajectory.
		strag_targ(Abeam,Zbeam,Zdepth,thetaBeam,phiBeam,Ebeam,thetaReact,phiReact,targRadL); 

		// Excitation energies
		// Select the excitation energies of the particles
		//Q = gsQvalue; 
		ExEject = 0.0; 
		ExRecoil = 0.0; 

		//ran_Ex_recoil = frand(); 
		//RecoilState = 0; 

		// Not the most efficient method in the world...
		if(ADist_type[RecoilState] == 0){// Isotropic, even weighting
			RecoilState = short(frand()*NRecoilStates);
			ExRecoil = ExRecoilStates[RecoilState];
		}
		else if(ADist_type[RecoilState] == 1){ // Angular dist weighted
			ran_Ex_recoil = ran_Ex_recoil*totXsect[NRecoilStates]; 
			for(i = 0; i < NRecoilStates; i++){
				if(ran_Ex_recoil > totXsect[i-1] && ran_Ex_recoil <= totXsect[i]){
					ExRecoil = ExRecoilStates[i]; 
					RecoilState = i; 
				} 
			} 

			for(j = 0; j < NAngPoints[RecoilState]; j++){ // set the angular distribution
				DistAng_temp[j] = DistAng[RecoilState][j]; 
				DistInt_temp[j] = DistInt[RecoilState][j]; 
			} 
		} 

		// the 2 body kinematics routine to generate the ejectile and recoil
		kindeux(Ereact,thetaReact,phiReact,gsQvalue,ExEject,ExRecoil,Abeam,Atarg,Arecoil,Aeject,ADist_type[RecoilState],NAngPoints[RecoilState],
			maxInt[RecoilState],DistAng_temp,DistInt_temp,thetaRecoil,phiRecoil,thetaEject,phiEject,Erecoil,Eeject);

		//file22 << thetaReact << "\t" << phiReact << "\n";
		//file22 << thetaEject << "\t" << Eeject << "\n";
		//file22 << phiEject << "\t" << Eeject << "\n";
		//std::cout << thetaEject << " " << phiEject << std::endl;

		// End of reaction calculation
	
		//---------------------------------------------------------------------------
		// Charged particle processing
		//---------------------------------------------------------------------------
		// Process the ejectile
		stuck = false; 
		if(Det_eject){
			// the charged reaction product routine to simulate the particle leaving the
			// target, check if it is detected and, if so, simulate the detection	
			Eject_hitGeo = 0;
			product_proc(Eeject,thetaEject,phiEject,Aeject,Zeject,beamspot,thickness,targ_depth,targ_angle,eject_E,eject_Erange,targRadL,
				     NDet_cyl,Nstrips_cyl,DetPhiMin,DetPhiMax,DetThetaMin,DetThetaMax,detLength_cyl,detWidth_cyl,DetRad_cyl,DetZ_cyl,
				     DetPhi_cyl,Ndet_plan,Nstrips_plan,DetZ_plan,DetXMax,DetXMin,DetYMax,DetYMin,Ndet_ann,Nstrips_ann,DetPhiMin_ann,
				     DetPhiMax_ann,DetThetaMin_ann,DetThetaMax_ann,eject_E_det,eject_Erange_det,conv_det,DetRadL,dEthick_cyl,
				     dEthick_plan,dEthick_ann,Ethick_cyl,Ethick_plan,Ethick_ann,Eres_dE_cyl,Pres_dE_cyl,Eres_E_cyl,Pres_E_cyl,
				     Pres_dE_en_cyl,Pres_E_en_cyl,Eres_dE_plan,Pres_dE_plan,Eres_E_plan,Pres_E_plan,Pres_dE_en_plan,Pres_E_en_plan,
				     xPlane,yPlane,detWidth_plan,detLength_plan,Eres_dE_ann,Eres_E_ann,DetZ_ann,DetInner_ann,DetOuter_ann,detPhi_ann,
				     EjectHit,Eject_hitGeo,DetHitCyl_eject,StripHitCyl_eject,DetHitPlan_eject,StripHitPlan_eject,DetHitAnn_eject,
				     StripHitAnn_eject,Eemerge_eject,thetaEjectdE,PhiEjectdE,det_dE_eject,det_E_eject,ThetaEjectE,PhiEjectE,
				     dE_final_eject,E_final_eject,dE_final_Theta_eject,dE_final_Phi_eject,E_final_Theta_eject,E_final_Phi_eject,stuck);
			
			// Set the ejecitle direction to that after straggling in the target
			// (NB - thetaEjectdE does not contain dE detector position resolution effects)
			if(!stuck){
				thetaEject = thetaEjectdE; 
				phiEject = PhiEjectdE;
			} 
		} // Det_eject

		//---------------------------------------------------------------------------
		// Process the RECOIL
		if(!stuck && Det_recoil){
			// the charged reaction product routine to simulate the particle leaving the
			// target, check if it is detected and, if so, simulate the detection
			Recoil_hitGeo = 0; 
			product_proc(Erecoil,thetaRecoil,phiRecoil,Arecoil,Zrecoil,beamspot,thickness,targ_depth,targ_angle,recoil_E,recoil_Erange,
				     targRadL,NDet_cyl,Nstrips_cyl,DetPhiMin,DetPhiMax,DetThetaMin,DetThetaMax,detLength_cyl,detWidth_cyl,DetRad_cyl,
				     DetZ_cyl,DetPhi_cyl,Ndet_plan,Nstrips_plan,DetZ_plan,DetXMax,DetXMin,DetYMax,DetYMin,Ndet_ann,Nstrips_ann,
				     DetPhiMin_ann,DetPhiMax_ann,DetThetaMin_ann,DetThetaMax_ann,recoil_E_det,recoil_Erange_det,conv_det,DetRadL,
				     dEthick_cyl,dEthick_plan,dEthick_ann,Ethick_cyl,Ethick_plan,Ethick_ann,Eres_dE_cyl,Pres_dE_cyl,Eres_E_cyl,
				     Pres_E_cyl,Pres_dE_en_cyl,Pres_E_en_cyl,Eres_dE_plan,Pres_dE_plan,Eres_E_plan,Pres_E_plan,Pres_dE_en_plan,
				     Pres_E_en_plan,xPlane,yPlane,detWidth_plan,detLength_plan,Eres_dE_ann,Eres_E_ann,DetZ_ann,DetInner_ann,
				     DetOuter_ann,detPhi_ann,RecoilHit,Recoil_hitGeo,DetHitCyl_recoil,StripHitCyl_recoil,DetHitPlan_recoil,
				     StripHitPlan_recoil,DetHitAnn_recoil,StripHitAnn_recoil,Eemerge_recoil,thetaRecoildE,PhiRecoildE,
				     det_dE_recoil,det_E_recoil,thetaRecoilE,PhiRecoilE,dE_final_recoil,E_final_recoil,dE_final_Theta_recoil,
				     dE_final_Phi_recoil,E_final_Theta_recoil,E_final_Phi_recoil,stuck);
				     
			// Set the recoil direction to that after straggling in the target
			// (NB - thetaRecoildE does not contain dE detector position resolution effects)
			if(!stuck){
				thetaRecoil = thetaRecoildE; 
				phiRecoil = PhiRecoildE; 
			}
		} // Det_recoil

		//------------------------------------------------------------------
		// OUTPUT FILES
		//------------------------------------------------------------------

		if(!stuck){
			// E-theta plot for all events
			// XXXXXX Put in logical to decide whether to do this, as it probably slows the program down XXXXXXX
			if(write_all_eject){ file11 << thetaEject*rad2deg << "\t" << Eeject << "\n"; } //(Eemerge_eject)

			//       if(write_all_recoil)
			//     &  write(41,*)(ThetaRecoil*rad2deg),(Erecoil) !(Eemerge_recoil)
			//	write(11,*)(E_final_Theta*rad2deg),(E_final_eject) ! Resolution
			//      write(13,*)(PhiEject*rad2deg),(ThetaEject*rad2deg) !Distribution of all events
			testArray[short(thetaEject*rad2deg)][0] += 1; 

			// Determine what's written out
			if(coincidence){// For coincidence hits
				if(EjectHit == 1 && RecoilHit == 1){
					write_eject = true; 
					write_recoil = true; 
				} 
				else{ 
					write_eject = false ; 
					write_recoil = false ; 
				} 
			} 
			else{ // For singles hits
				if(EjectHit == 1){ write_eject = true; } 
				else{ write_eject = false; } 
				if(RecoilHit == 1){ write_recoil = true; } 
				else{ write_recoil = false; } 
			} 

			//------------------------------------------------------------------
			// For detected EJECTILES
			if(write_eject){
				// E-theta plot for detected events,Unsmeared
				file12 << (thetaEject*rad2deg) << (Eemerge_eject) << "\n"; 

				// Coverage in theta-phi space, Unsmeared
				file13 << (phiEject*rad2deg) << (thetaEject*rad2deg) << "\n"; 

				// Etot-thetaE plot for detected events, smeared
				if(E_final_eject > 0.0){ file14 << E_final_Theta_eject*rad2deg << "\t" << E_final_eject + dE_final_eject << "\n"; }

				// Etot-thetadE plot for detected events, smeared
				file15 << (dE_final_Theta_eject*rad2deg) << (E_final_eject + dE_final_eject) << "\n"; 

				// dE-thetaE plot for detected events, smeared
				if(E_final_eject > 0.0){ file16 << E_final_Theta_eject*rad2deg << "\t" << dE_final_eject << "\n"; }

				// dE-thetadE plot for detected events, smeared
				file17 << (dE_final_Theta_eject*rad2deg) << (dE_final_eject) << "\n"; 

				// E_residual-thetaE plot for detected events, smeared
				if(E_final_eject > 0.0){ file18 << E_final_Theta_eject*rad2deg << "\t" << E_final_eject << "\n"; }

				// E_residual-thetadE plot for detected events, smeared
				if(E_final_eject > 0.0){ file19 << dE_final_Theta_eject*rad2deg << "\t" << E_final_eject << "\n"; }
			
				if(thetaEject <= pi){ file20 << (E_final_eject) << (dE_final_eject) << "\n"; } // PI plot - forward angles
				else{ file21 << (E_final_eject) << (dE_final_eject) << "\n"; } // PI plot - backward angles

			} // write_eject

			//------------------------------------------------------------------
			// For detected RECOILS
			//       if(write_recoil)then
			// E-theta plot for detected events,Unsmeared
			//	  write(42,*)(ThetaRecoil*rad2deg),(Eemerge_recoil)
			//
			// Coverage in theta-phi space, Unsmeared
			//	 write(43,*)(PhiRecoil*rad2deg),(ThetaRecoil*rad2deg)
			//
			// Etot-thetaE plot for detected events, smeared
			//	  if(E_final_recoil.gt.0.0)
			//     &       write(44,*)(E_final_theta_recoil*rad2deg),
			//     &       (E_final_recoil + dE_final_recoil)
			// Etot-thetadE plot for detected events, smeared
			//	    write(45,*)(dE_final_Theta_recoil*rad2deg),
			//     &       (E_final_recoil + dE_final_recoil)
			//
			// dE-thetaE plot for detected events, smeared
			//	  if(E_final_recoil.gt.0.0)
			//     &       write(46,*)(E_final_theta_recoil*rad2deg),
			//     &       (dE_final_recoil)
			//
			// dE-thetadE plot for detected events, smeared
			//	 write(47,*)(dE_final_Theta_recoil*rad2deg),
			//     &   (dE_final_recoil)
			//
			// E_residual-thetaE plot for detected events, smeared
			//	 if(E_final_recoil.gt.0.0) write(48,*)
			//     &    (E_final_theta_recoil*rad2deg),(E_final_recoil)
			//
			// E_residual-thetadE plot for detected events, smeared
			//	 if(E_final_recoil.gt.0.0) write(49,*)
			//     &  (dE_final_Theta_recoil*rad2deg),(E_final_recoil)
			//
			// PI plot - forward angles
			//	  if(ThetaRecoil.le.pi)then
			//	   write(22,*) (E_final_recoil),(dE_final_recoil)
			// PI plot - backward angles
			//	  else
			//	    write(23,*) (E_final_recoil),(dE_final_recoil)
			//	  endif
			//	endif ! write_recoil
			//
			// Storage for output
			// Store the energy-angle plot as a 2D spectrum for EJECTILES
			//      eject_Et(nshort(thetaEject/pi*360.0),nshort(Eemerge_eject*10.0)) = 
			//     & eject_Et(nshort(thetaEject/pi*360.0),nshort(Eemerge_eject*10.0))+1
			//
			// Store the energy-angle plot as a 2D spectrum for RECOILS
			//      recoil_Et(nshort(thetaRecoil/pi*360.0),nshort(Eemerge_recoil*0.5)) = 
			//     & recoil_Et(nshort(thetaRecoil/pi*360.0),nshort(Eemerge_recoil*0.5))+1
			//	std::cout << "'E_recoil = ', Erecoil
			//	std::cout << "'E_ejectile = ', Eeject
			//
			// smearing of the reaction energy & location throughout the target
			//
			// smearing of the angular spread of beam particles throughout the target
			// smearing of the ejectiles throughout the target.

			// ###########################################################################
			//          Ndetected = Ndetected + 1
		
			if(coincidence){
				if((Det_eject && EjectHit == 1) && (Det_recoil && RecoilHit == 1)){ 
					Ndetected = Ndetected + 1; 
					if(!flag){ flag = true; }
				}
			}
			else if((Det_eject && EjectHit == 1) || (Det_recoil && RecoilHit == 1)){ 
				Ndetected = Ndetected + 1; 
				if(!flag){ flag = true; }
			}
		}
		Nsimulated = Nsimulated + 1; 

	} // Main simulation loop
	// ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  = 

	std::ofstream file31("testArray.dat"); 
	for(i = 0; i < 360; i++){
		//        do j = 0,360
		//	 if(testArray(i,j).gt.0) write(21,*) i,j,testArray(i,j)
		if(testArray[i][0]>0){ file31 << i << "\t" << testArray[i][0] << "\n"; }
		//	enddo
	} 

	file11.close(); 
	file12.close(); 
	file13.close(); 
	file14.close(); 
	file15.close(); 
	file16.close(); 
	file17.close(); 
	file18.close(); 
	file19.close(); 
	file20.close(); 
	file21.close(); 
	file22.close(); // (CORY)
	//      close(22)
	//      close(23)
	//      close(41)
	//      close(42)
	//      close(43)
	//      close(44)
	//      close(45)
	//      close(46)
	//      close(47)
	//      close(48)
	//      close(49)
	//      close(57)!debug
	//      close(58)!debug

	file31.close(); 
 
	std::cout << "\n ------------- Simulation Complete --------------\n";
	std::cout << " The simulation took " << (float)(clock()-timer)/CLOCKS_PER_SEC << " seconds\n\n"; 
} 
