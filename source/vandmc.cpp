/** \file vandmc.cpp
 * \brief Ray tracing kinematics simulations program for arbitrary setups.
 *
 * This program is intended to be used as a way of performing kinematics 
 * simulations on arbitary combinations of detector geometries. It should 
 * provide a semi-transparent way to visualize what is happening in a low 
 * energy reaction. It is written in the spirit of VIKAR by Steve Pain.
 *
 * \author C. R. Thornsberry
 * \date Feb. 26th, 2016
 */
#include <fstream>
#include <iostream>
#include <time.h>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"

// SimpleScan
#include "CTerminal.h"

// VANDMC
#include "vandmc.hpp"

#define VERSION "1.33"

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

///////////////////////////////////////////////////////////////////////////////
// class vandmcParameter
///////////////////////////////////////////////////////////////////////////////

double vandmcParameter::ConvertToDouble(){
	return strtod(val.c_str(), NULL);
}

unsigned int vandmcParameter::ConvertToUlong(){
	return strtoul(val.c_str(), NULL, 10);
}

int vandmcParameter::ConvertToLong(){
	return strtol(val.c_str(), NULL, 10);
}

bool vandmcParameter::ConvertToBool(){
	return (bool)strtol(val.c_str(), NULL, 10);
}

///////////////////////////////////////////////////////////////////////////////
// class vandmcParameterReader
///////////////////////////////////////////////////////////////////////////////

bool vandmcParameterReader::Read(const char *fname, bool echo/*=false*/){
	std::ifstream ifile(fname);
	if(!ifile.good()) return false;

	std::string line;
	while(true){
		getline(ifile, line);
		if(ifile.eof()) break;
	
		// Check for blank lines or comment lines.
		if(line.empty() || line[0] == '#') continue;

		// Split the string into parts.
		std::vector<std::string> args;
		split_str(line, args, '\t');

		// Check for invalid number of arguments.
		if(args.size() < 2) continue;

		bool foundParamName;
		for(std::vector<vandmcParameter>::iterator iter = validParameters.begin(); iter != validParameters.end(); iter++){
			if(iter->Compare(args.at(0))){
				if(echo) std::cout << args.at(0) << "=" << args.at(1) << std::endl;
				parameters.push_back(vandmcParameter(args.at(0), "", args.at(1)));
				foundParamName = true;
				break;
			}
		}

		if(!foundParamName) std::cout << " vandmc: Encountered unrecognized parameter name (" << args.at(0) << ")\n";
	}
	ifile.close();	

	return true;
}

bool vandmcParameterReader::FindInList(const std::string &str){
	return (findParam(str) != NULL);
}

size_t vandmcParameterReader::FindAllOccurances(const std::string &str, std::vector<vandmcParameter*> &vec){
	vec.clear();
	for(std::vector<vandmcParameter>::iterator iter = parameters.begin(); iter != parameters.end(); iter++){
		if(iter->Compare(str) && !iter->Empty()) vec.push_back(&(*iter));
	}
	return vec.size();
}

bool vandmcParameterReader::FindDouble(const std::string &str, double &val){
	vandmcParameter *par = findParam(str);
	if(par == NULL) return false;
	val = par->ConvertToDouble();
	return true;
}

bool vandmcParameterReader::FindUlong(const std::string &str, unsigned int &val){
	vandmcParameter *par = findParam(str);
	if(par == NULL) return false;
	val = par->ConvertToUlong();
	return true;
}

bool vandmcParameterReader::FindLong(const std::string &str, int &val){
	vandmcParameter *par = findParam(str);
	if(par == NULL) return false;
	val = par->ConvertToLong();
	return true;
}

bool vandmcParameterReader::FindBool(const std::string &str, bool &val){
	vandmcParameter *par = findParam(str);
	if(par == NULL) return false;
	val = par->ConvertToBool();
	return true;
}

bool vandmcParameterReader::FindString(const std::string &str, std::string &val){
	vandmcParameter *par = findParam(str);
	if(par == NULL) return false;
	val = par->GetValue();
	return true;
}

vandmcParameter *vandmcParameterReader::findParam(const std::string &str){
	for(std::vector<vandmcParameter>::iterator iter = parameters.begin(); iter != parameters.end(); iter++){
		if(iter->Compare(str)){
			if(iter->Empty()) return NULL;
			return &(*iter);
		}
	}
	return NULL;
}

void vandmcParameterReader::initialize(){
	const std::vector<std::string> paramNames = {"VERSION",
	                                             "BEAM_Z",
	                                             "BEAM_A",
	                                             "BEAM_AMU",
	                                             "TARG_Z",
	                                             "TARG_A",
	                                             "TARG_AMU",
	                                             "TARG_MATERIAL",
	                                             "TARG_THICKNESS",
	                                             "TARG_ANGLE",
	                                             "RECOIL_Z",
	                                             "RECOIL_A",
	                                             "RECOIL_AMU",
	                                             "EJECT_Z",
	                                             "EJECT_A",
	                                             "EJECT_AMU",
	                                             "BEAM_E",
	                                             "BEAM_TYPE",
	                                             "BEAM_SPOT",
	                                             "BEAM_ANGLE",
	                                             "BEAM_RATE",
	                                             "BEAM_E_SPREAD",
	                                             "RECOIL_STATE",
	                                             "ANGULAR_DIST_MODE",
	                                             "ANGULAR_DIST",
	                                             "SMALL_EFFICIENCY",
	                                             "MED_EFFICIENCY",
	                                             "LARGE_EFFICIENCY",
	                                             "DETECTOR_FNAME",
	                                             "N_SIMULATED_PART",
	                                             "BACKGROUND_RATE",
	                                             "BACKGROUND_WINDOW",
	                                             "BACKGROUND_WAIT",
	                                             "REQUIRE_COINCIDENCE",
	                                             "WRITE_REACTION_INFO",
	                                             "SIMULATE_252CF"};

	for(std::vector<std::string>::const_iterator iter = paramNames.begin(); iter != paramNames.end(); iter++){
		validParameters.push_back(vandmcParameter(*iter));
	}
}

///////////////////////////////////////////////////////////////////////////////
// class vandmc
///////////////////////////////////////////////////////////////////////////////

void vandmc::initialize(){
	// Seed randomizer
	srand(time(NULL));

	num_materials = 0;
	
	// Physics Variables
	NRecoilStates = 0;
	ExRecoilStates = NULL;
	totXsect = NULL; 
	gsQvalue = 0.0;

	// Energy variables
	Ebeam = 0.0, Ebeam0 = 0.0;
	ErecoilMod = 0.0;
	EejectMod = 0.0;
	Egamma = 0.0;
		
	// Beam variables
	beamspot = 0.0; // Beamspot diameter (m) (on the surface of the target)
	beamEspread = 0.0; // Beam energy spread (MeV)
	beamAngdiv = 0.0; // Beam angular divergence (radians)

	timeRes = 2E-9; // Pixie-16 time resolution (s)
	BeamRate = 0.0; // Beam rate (1/s)

	// Background variables
	backgroundRate = 0;
	backgroundWait = 0;
	detWindow = 0.0;
	bgPerDetection = false;
	
	// Input/output variables
	NgoodDetections = 0; // Number of wanted particles detected in ejectile detectors.
	Ndetected = 0; // Total number of particles detected in ejectile detectors.
	Nwanted = 0; // Number of desired ejectile detections.
	Nsimulated = 0; // Total number of simulated particles.
	NdetHit = 0; // Total number of particles which collided with a detector.
	Nreactions = 0; // Total number of particles which react with the target.
	NrecoilHits = 0;
	NejectileHits = 0;
	NgammaHits = 0;
	NvetoEvents = 0;
	Ndet = 0; // Total number of detectors
	NdetRecoil = 0; // Total number of recoil detectors
	NdetEject = 0; // Total number of ejectile detectors
	NdetGamma = 0; // Total number of gamma detectors
	NdetVeto = 0; // Total number of particle vetos
	BeamType = 0; // The type of beam to simulate (0=gaussian, 1=cylindrical, 2=halo)

	// Default options
	InverseKinematics = true;
	InCoincidence = true;
	WriteReaction = false;
	NeutronSource = false;
	PerfectDet = true;
	SupplyRates = false;
	BeamFocus = false;
	DoRutherford = false;
	echoMode = false;
	printParams = false;
	ADists = 0;
	
	// Detector options
	have_recoil_det = false;
	have_ejectile_det = false;
	have_gamma_det = false;
	have_veto_det = false;

	// Output filename string
	output_filename = "vandmc.root";

	handler.add(optionExt("input", required_argument, NULL, 'i', "<filename>", "Specify an input configuration file."));
	handler.add(optionExt("output", required_argument, NULL, 'o', "<filename>", "Specify the name of the output file."));
	handler.add(optionExt("detector", required_argument, NULL, 'd', "<filename>", "Specify the name of the detector file."));
	handler.add(optionExt("echo", no_argument, NULL, 'e', "", "Echo values read from the config file."));
	handler.add(optionExt("print", no_argument, NULL, 0x0, "", "Print simulation parameters."));
}

void vandmc::titleCard(){
	std::cout << "####    ####      ##      ####    #### ########     ####    ####  #########   \n";
	std::cout << " ##      ##       ##       ##      ##   ##   ###     ##      ##  ###     ###  \n";
	std::cout << "  ##    ##       ####      ###     ##   ##     ###   ###    ###  ##        ## \n";
	std::cout << "  ##    ##       ####      ####    ##   ##       ##  ####  ####  ##           \n";
	std::cout << "   ##  ##       ##  ##     ## ##   ##   ##       ##  ## #### ##  ##           \n";
	std::cout << "   ##  ##       ##  ##     ##  ##  ##   ##       ##  ##  ##  ##  ##           \n";
	std::cout << "    ####       ########    ##   #####   ##       ##  ##      ##  ##           \n";
	std::cout << "    ####       ##    ##    ##    ####   ##     ###   ##      ##  ##        ## \n";
	std::cout << "     ##       ##      ##   ##      ##   ##   ###     ##      ##  ###     ###  \n";
	std::cout << "     ##      ####    #### ####    #### ########     ####    ####  #########   \n";

	std::cout << "\n VANDMC v " << VERSION << "\n"; 
	std::cout << " ==  ==  ==  ==  == \n\n"; 
	
	std::cout << " Welcome to VANDMC, the Versatile Array of Neutron Detectors Monte Carlo program.\n";
	std::cout << "  Loading the input configuration file...\n";
	
	std::cout << "\n ==  ==  ==  ==  == \n\n";
	sleep(1);
}

bool vandmc::readConfig(const char *fname){
	if(!reader.Read(fname, echoMode)) return false;

	if(echoMode) std::cout << std::endl;

	std::cout << " Reading from file " << fname << std::endl;

	std::string str;
	double dval;
	
	std::string versionString;
	if(reader.FindString("VERSION", versionString))
		std::cout << "  Input file version: " << versionString << std::endl;

	if(reader.FindDouble("BEAM_Z", dval))
		beam_part.SetZ(dval);
	if(reader.FindDouble("BEAM_A", dval))
		beam_part.SetA(dval);	
	if(reader.FindDouble("BEAM_AMU", dval))
		beam_part.SetMass(dval);	

	if(reader.FindDouble("TARG_Z", dval))
		targ.SetZ(dval);
	if(reader.FindDouble("TARG_A", dval))
		targ.SetA(dval);	
	if(reader.FindDouble("TARG_AMU", dval))
		targ.SetMass(dval);		

	if(reader.FindDouble("RECOIL_Z", dval))
		recoil_part.SetZ(dval);
	if(reader.FindDouble("RECOIL_A", dval))
		recoil_part.SetA(dval);	
	if(reader.FindDouble("RECOIL_AMU", dval))
		recoil_part.SetMass(dval);	

	if(reader.FindDouble("EJECT_Z", dval))
		eject_part.SetZ(dval);
	if(reader.FindDouble("EJECT_A", dval))
		eject_part.SetA(dval);	
	if(reader.FindDouble("EJECT_AMU", dval))
		eject_part.SetMass(dval);	

	// Set the reaction Q-value.
	gsQvalue = (beam_part.GetMass()+targ.GetMass())-(recoil_part.GetMass()+eject_part.GetMass()); 
	
	// Set inverse or normal kinematics.
	InverseKinematics = (beam_part.GetA() > targ.GetA() ? true : false);

	reader.FindDouble("BEAM_E", Ebeam0); // In MeV
	reader.FindUlong("BEAM_TYPE", BeamType);
	if(BeamType > 2){ 
		std::cout << " FATAL ERROR! Invalid beam type selection (" << BeamType << ").\n";
		return false;
	}
	if(reader.FindDouble("BEAM_SPOT", dval))
		beamspot = dval/1000.0; // in meters
	if(reader.FindDouble("BEAM_ANGLE", dval))
		beamAngdiv = abs(dval*deg2rad); // in radians
	reader.FindDouble("BEAM_E_SPREAD", beamEspread); // In MeV
	
	std::vector<vandmcParameter*> tempParams;
	NRecoilStates = reader.FindAllOccurances("RECOIL_STATE", tempParams);
	ExRecoilStates = new double[NRecoilStates+1]; ExRecoilStates[0] = 0.0;
	totXsect = new double[NRecoilStates+1]; totXsect[0] = 0.0;
	
	size_t state = 1;
	for(std::vector<vandmcParameter*>::iterator iter = tempParams.begin(); iter != tempParams.end(); iter++){
		ExRecoilStates[state] = (*iter)->ConvertToDouble();
		totXsect[state++] = 0.0;
	}
	
	// User provided angular distributions.
	reader.FindUlong("ANGULAR_DIST_MODE", ADists);
	if(ADists == 1 || ADists == 2){ // Read the filenames.
		unsigned int NAngularDists = reader.FindAllOccurances("ANGULAR_DIST", tempParams);
		if(NAngularDists != NRecoilStates){
			std::cout << " FATAL ERROR! Too few occurances of \"ANGULAR_DIST\" in config file (" << NAngularDists << " < " << NRecoilStates << ").\n";
			return false;
		}
		
		if(ADists == 1){ // Read in filename for each state's angular distribution.
			for(std::vector<vandmcParameter*>::iterator iter = tempParams.begin(); iter != tempParams.end(); iter++){
				if(!DoRutherford && (*iter)->GetValue() != "RUTHERFORD")
					AngDist_fname.push_back((*iter)->GetValue());
				else
					DoRutherford = true;
			}

			// Supply beam rate information.
			if(reader.FindDouble("BEAM_RATE", dval)){
				SupplyRates = true;
				BeamRate = dval; // in pps
			}
		}
		else{ // Supply the rates relative to the ground state. 
			// Example: a value of 1.0 will generate 1 particle for each ground state particle. A
			// value of 2.0 will generate 2 particles for each ground state and so on.
			for(std::vector<vandmcParameter*>::iterator iter = tempParams.begin(); iter != tempParams.end(); iter++){
				AngDist_fname.push_back((*iter)->GetValue());
			}
		}
	}

	// Target material name
	reader.FindString("TARG_MATERIAL", targ_mat_name);
		
	// Target thickness
	if(reader.FindDouble("TARG_THICKNESS", dval)){
		if(targ_mat_name != "NONE") // Normal material.
			targ.SetThickness(dval);
		else // Target material energy loss disabled.
			targ.SetRealThickness(dval);
	}
		
	// Target angle wrt beam axis
	if(reader.FindDouble("TARG_ANGLE", dval))
		targ.SetAngle(dval*deg2rad);

	// Detector efficiencies
	if(reader.FindString("SMALL_EFFICIENCY", str)){
		bar_eff.ReadSmall(str.c_str());
		PerfectDet = false;
	}
	if(reader.FindString("MED_EFFICIENCY", str)){
		bar_eff.ReadMedium(str.c_str());
		PerfectDet = false;
	}
	if(reader.FindString("LARGE_EFFICIENCY", str)){
		bar_eff.ReadLarge(str.c_str());
		PerfectDet = false;
	}

	// Detector setup filename
	reader.FindString("DETECTOR_FNAME", detector_filename);

	// Maximum number of detected particles
	reader.FindUlong("N_SIMULATED_PART", Nwanted);
	
	// Background rate (per detection event)
	reader.FindUlong("BACKGROUND_RATE", backgroundRate);
		
	if(backgroundRate > 0){ 
		// Get the detection ToF window
		reader.FindDouble("BACKGROUND_WINDOW", detWindow); // in ns
			
		backgroundWait = backgroundRate;

		reader.FindBool("BACKGROUND_PER_RECOIL", bgPerDetection);
	}

	reader.FindBool("REQUIRE_COINCIDENCE", InCoincidence);
	reader.FindBool("WRITE_REACTION_INFO", WriteReaction);
	reader.FindBool("SIMULATE_252CF", NeutronSource);

	return true;
}

bool vandmc::setup(int argc, char *argv[]){
	if(!handler.setup(argc, argv)){
		return false;
	}

	// Set input filename
	if(!handler.getOption(0)->active){
		std::cout << " FATAL ERROR! No input filename specified!\n";
		return false;
	}
	input_filename = handler.getOption(0)->argument;

	// Set output filename
	if(handler.getOption(1)->active){
		output_filename = handler.getOption(1)->argument;
	}

	// Set detector filename
	if(handler.getOption(2)->active){
		detector_filename = handler.getOption(2)->argument;
	}

	// Set config file read echo mode
	if(handler.getOption(3)->active){
		echoMode = true;
	}

	// Set parameter print mode
	if(handler.getOption(4)->active){
		printParams = true;
	}	
	

	return true;
}

void vandmc::print(){
	std::cout << "  Beam-Z: " << beam_part.GetZ() << std::endl;
	std::cout << "  Beam-A: " << beam_part.GetA() << std::endl;
	std::cout << "  Beam Mass: " << beam_part.GetMassAMU() << " amu\n";
	std::cout << "  Target-Z: " << targ.GetZ() << std::endl;
	std::cout << "  Target-A: " << targ.GetA() << std::endl;
	std::cout << "  Target Mass: " << targ.GetMassAMU() << " amu\n";
	std::cout << "  Recoil-Z: " << recoil_part.GetZ() << std::endl;
	std::cout << "  Recoil-A: " << recoil_part.GetA() << std::endl;
	std::cout << "  Recoil Mass: " << recoil_part.GetMassAMU() << " amu\n";
	std::cout << "  Ejectile-Z: " << eject_part.GetZ() << std::endl;
	std::cout << "  Ejectile-A: " << eject_part.GetA() << std::endl;
	std::cout << "  Ejectile Mass: " << eject_part.GetMassAMU() << " amu\n";
	std::cout << "  G.S. Q-Value: " << gsQvalue << " MeV\n";
	std::cout << "  Inverse Kinematics: " << (InverseKinematics ? "YES" : "NO") << "\n";
	std::cout << "  Beam Energy: " << Ebeam0 << " MeV\n";
	if(BeamType == 0)
		std::cout << "  Beam Spot FWHM: " << beamspot*1000 << " mm\n"; 
	else if(BeamType == 1)
		std::cout << "  Beam Spot Diameter: " << beamspot*1000 << " mm\n"; 
	else if(BeamType == 2)
		std::cout << "  Beam Halo Diameter: " << beamspot*1000 << " mm\n"; 
	std::cout << "  Beam Angular Divergence: " << beamAngdiv*rad2deg << " degrees\n";
	std::cout << "  Beam Energy Spread: " << beamEspread << " MeV\n";
	std::cout << "  No. Excited States: " << NRecoilStates-1 << std::endl;
	std::cout << "   Recoil Ground State: 0.0 MeV\n";
	for(unsigned int i = 1; i < NRecoilStates; i++)
		std::cout << "   Recoil Excited State " << i << ": " << ExRecoilStates[i] << " MeV\n";
	std::cout << "  Supply Angular Distributions: " << (ADists == 1 || ADists == 2 ? "YES" : "NO") << "\n";
	if(ADists == 1){
		if(!DoRutherford){
			for(unsigned int i = 0; i < NRecoilStates; i++){
				if(i == 0)
					std::cout << "   Distribution for ground state: " << AngDist_fname.at(i) << std::endl;
				else
					std::cout << "   Distribution for state " << i+1 << ": " << AngDist_fname.at(i) << std::endl;
			}
		}
		else{
			for(unsigned int i = 0; i < NRecoilStates; i++){
				if(i == 0)
					std::cout << "   Distribution for ground state: RUTHERFORD\n";
				else
					std::cout << "   Distribution for state " << i+1 << ": RUTHERFORD\n";
			}
		}
		std::cout << "   Beam Rate: " << BeamRate << " pps\n";
	}
	else if(ADists == 2){
		for(unsigned int i = 0; i < NRecoilStates; i++){
			if(i == 0)
				std::cout << "   Production rate ground state: " << AngDist_fname.at(i) << " per event.\n";
			else 
				std::cout << "   Production rate for state " << i+1 << ": " << AngDist_fname.at(i) << " per event.\n";
		}
	}
	std::cout << "  Target Material: " << targ_mat_name << std::endl;
	if(targ_mat_name != "NONE")
		std::cout << "  Target Thickness: " << targ.GetThickness() << " mg/cm^2\n";	
	else
		std::cout << "  Target Thickness: " << targ.GetRealThickness() << " m\n";	
	std::cout << "  Target Angle: " << targ.GetAngle()*rad2deg << " degrees\n";
	std::cout << "  Perfect Detectors: " << (PerfectDet ? "YES" : "NO") << "\n";
	if(bar_eff.GetNsmall() > 0)
		std::cout << "   Found " << bar_eff.GetNsmall() << " small bar efficiency data points.\n";
	if(bar_eff.GetNmedium() > 0)
		std::cout << "   Found " << bar_eff.GetNmedium() << " medium bar efficiency data points.\n";
	if(bar_eff.GetNlarge() > 0)
		std::cout << "   Found " << bar_eff.GetNlarge() << " large bar efficiency data points.\n";
	std::cout << "  Detector Setup Filename: " << detector_filename << std::endl;
	std::cout << "  Desired Detections: " << Nwanted << std::endl; 
	if(backgroundRate > 0){
		if(bgPerDetection)
			std::cout << "  Background Rate: " << backgroundRate << " events per detection\n";
		else
			std::cout << "  Background Rate: " << backgroundRate << " events per recoil\n";
		std::cout << "  Background Time Width: " << detWindow << " ns\n";
	}
	else
		std::cout << "  Background Rate: NONE\n";
	std::cout << "  Require Particle Coincidence: " << (InCoincidence ? "YES" : "NO") << std::endl;
	std::cout << "  Write Reaction Info: " << (WriteReaction ? "YES" : "NO") << std::endl;
	std::cout << "  Simulate 252Cf source: " << (NeutronSource ? "YES" : "NO") << std::endl;
}

bool vandmc::Execute(int argc, char *argv[]){ 
	// Set all variables to default values.
	initialize();

	// Handle all command line arguments
	if(!setup(argc, argv)) return false;

	// Display the vandmc welcome screen.
	titleCard();
	
	// Read the input config file.
	if(!readConfig(input_filename.c_str())){
		std::cout << " FATAL ERROR! Failed to read configuration file \"" << input_filename << "\"!\n";
		return false;
	}
	
	if(printParams) print();
	
	// Check for mass conservation
	if(beam_part.GetA()+targ.GetA() != recoil_part.GetA()+eject_part.GetA()){
		std::cout << "\n FATAL ERROR! Mass number is NOT conserved! Aborting\n";
		return false;
	}
	
	// Check for no detector file defined.
	if(detector_filename.empty()){
		std::cout << "\n FATAL ERROR! Detector file not specified.\n";
		return false;
	}
		
	std::cout << "\n ==  ==  ==  ==  == \n\n";

	// Make sure the input variables are correct
	if(!Prompt(" Are the above settings correct?")){
		std::cout << "  ABORTING...\n";
		return false;
	}

	std::cout << "\n Initializing main simulation Kindeux object...\n";

	// Initialize kinematics object
	kind.Initialize(beam_part.GetA(), targ.GetA(), recoil_part.GetA(), eject_part.GetA(), gsQvalue, NRecoilStates, ExRecoilStates);

	// Set the simulated 252Cf source.
	if(NeutronSource) kind.ToggleNeutronSource();

	// Read the detector setup file
	std::cout << " Reading in NewVANDMC detector setup file...\n";
	Ndet = ReadDetFile(detector_filename.c_str(), vandle_bars);
	if(Ndet < 0){ // Failed to load setup file
		std::cout << " FATAL ERROR! failed to load detector setup file!\n";
		return false; 
	}
	else if(Ndet == 0){ // Check there's at least 1 detector!
		std::cout << " FATAL ERROR! Found no detectors in the detector setup file!\n"; 
		return false;
	}
	
	std::vector<std::string> needed_materials;
	for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
		if(!IsInVector((*iter)->GetMaterialName(), needed_materials)){
			needed_materials.push_back((*iter)->GetMaterialName());
		}
		
		if((*iter)->IsEjectileDet()){ NdetEject++; }
		if((*iter)->IsRecoilDet()){ NdetRecoil++; }
		if((*iter)->IsGammaDet()){ NdetGamma++; }
		if((*iter)->IsVeto()){ NdetVeto++; }
	}

	if(NdetRecoil > 0){ have_recoil_det = true; }
	if(NdetEject > 0){ have_ejectile_det = true; }
	if(NdetGamma > 0){ have_gamma_det = true; }
	if(NdetVeto > 0){ have_veto_det = true; }

	if(InCoincidence && (NdetEject == 0 || NdetGamma == 0)){
		std::cout << " FATAL ERROR! Requiring particle coincidence but found no ejectile or gamma detectors in the detector setup file!\n"; 
		return false;
	}

	if(InCoincidence && NdetRecoil == 0){
		std::cout << " FATAL ERROR! Requiring particle coincidence but found no recoil detectors in the detector setup file!\n"; 
		return false;
	}

	// Report on how many detectors were read in
	std::cout << " Found the following detector types in detector setup file " << detector_filename << std::endl;
	std::cout << "  Ejectile: " << NdetEject << std::endl; 
	std::cout << "  Recoil:   " << NdetRecoil << std::endl;
	std::cout << "  Gamma:    " << NdetGamma << std::endl;
	std::cout << "  Veto:     " << NdetVeto << std::endl;

	if(!have_recoil_det && !have_ejectile_det && !have_gamma_det){
		std::cout << " FATAL ERROR! Found no valid detectors in detector setup file!\n";
		return false;
	}

	std::cout << "\n Setting up VANDMC materials...\n";

	// Natural Gold
	materials.push_back(Material("Au197", 19.311, 79, 196.96657, 1));
	// BC408 plastic scintillator (polyvinyltoluene (C9H10 118.18 g/mol) base)
	materials.push_back(Material("BC408", 1.032, 6, 12.0107, 9, 1, 1.00794, 10));
	// Deuterated polyethylene
	materials.push_back(Material("C2D4", 1.06300, 6, 12.0107, 2, 1, 2.01588, 4));
	// Polyethylene
	materials.push_back(Material("C2H4", 0.95, 6, 12.0107, 2, 1, 1.00794, 4));
	// Polystyrene
	materials.push_back(Material("C8H8", 1.05, 6, 12.0107, 8, 1, 1.00794, 8));
	// Natural Carbon
	materials.push_back(Material("C12", 2.2670, 6, 12.0107, 1));
	// Deuterated polyethylene
	materials.push_back(Material("CD2", 1.06300, 6, 12.0107, 1, 1, 2.01588, 2));	
	// Silicon
	materials.push_back(Material("Si28", 2.3212, 14, 28.0855, 1));

	num_materials = materials.size();
	std::cout << " Successfully setup " << num_materials << " materials\n";

	bool use_target_eloss = true;
	targ_mat_id = 0;
	
	if(targ_mat_name == "NONE"){ // Use no target.
		std::cout << "  Target Material: DISABLED\n";
		use_target_eloss = false;
	}
	else{ // Use a pre-loaded material target.
		bool found_valid_material = false;
		std::cout << "  Target Material: " << targ_mat_name;
		for(targ_mat_id = 0; targ_mat_id < num_materials; targ_mat_id++){
			if(targ_mat_name == materials[targ_mat_id].GetName()){
				found_valid_material = true;
				break;
			}
		}
		if(!found_valid_material){ 
			std::cout << " (not found)"; 
			use_target_eloss = false;			
		}
		std::cout << std::endl;
	}
	
	if(use_target_eloss){
		targ.SetDensity(materials[targ_mat_id].GetDensity());
		targ.SetRadLength(materials[targ_mat_id].GetRadLength());
		std::cout << "  Target Radiation Length: " << targ.GetRadLength() << " mg/cm^2\n\n";
	
		// Calculate the stopping power table for the reactino particles in the target
		if(beam_part.GetZ() > 0){ // The beam is a charged particle (not a neutron)
			std::cout << " Calculating range table for beam in " << materials[targ_mat_id].GetName() << "...";
			beam_targ.Init(1000, 0.1, (Ebeam0+2*beamEspread), beam_part.GetZ(), beam_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
		if(eject_part.GetZ() > 0){
			std::cout << " Calculating range table for ejectile in " << materials[targ_mat_id].GetName() << "...";
			eject_targ.Init(1000, 0.1, (Ebeam0+2*beamEspread), eject_part.GetZ(), eject_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
		if(recoil_part.GetZ() > 0){
			std::cout << " Calculating range table for recoil in " << materials[targ_mat_id].GetName() << "...";
			recoil_targ.Init(1000, 0.1, (Ebeam0+2*beamEspread), recoil_part.GetZ(), recoil_part.GetA()/mev2amu, &materials[targ_mat_id]);
			std::cout << " Done!\n";
		}
	}

	// Set the molar mass of the target.
	targ.SetMolarMass(materials[targ_mat_id].GetMolarMass());

	// Calculate the stopping power table for the ejectiles in the materials
	if(eject_part.GetZ() > 0){ // The ejectile is a charged particle (not a neutron)
		eject_tables.assign(num_materials, RangeTable());
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating ejectile range table for " << materials[i].GetName() << "...";
			eject_tables[i].Init(1000, eject_part.GetKEfromV(0.02*c), (Ebeam0+2*beamEspread), eject_part.GetZ(), eject_part.GetA()/mev2amu, &materials[i]);
			std::cout << " Done!\n";
		}
	}

	// Calculate the stopping power table for the recoils in the materials
	if(recoil_part.GetZ() > 0){ // The recoil is a charged particle (not a neutron)
		recoil_tables.assign(num_materials, RangeTable());
		for(unsigned int i = 0; i < num_materials; i++){
			if(!IsInVector(materials[i].GetName(), needed_materials)){ continue; }
			std::cout << " Calculating recoil range table for " << materials[i].GetName() << "...";
			recoil_tables[i].Init(1000, recoil_part.GetKEfromV(0.02*c), (Ebeam0+2*beamEspread), recoil_part.GetZ(), recoil_part.GetA()/mev2amu, &materials[i]);
			std::cout << " Done!\n";
		}
	}
	
	// Calculate the beam focal point (if it exists)
	lab_beam_focus = Vector3(0.0, 0.0, 0.0);
	if(beamAngdiv >= 0.000174532925199){
		lab_beam_focus.axis[2] = -beamspot/(2.0*std::tan(beamAngdiv));
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
			if(kind.SetDist(AngDist_fname, BeamRate, &targ)){
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
		return false;
	}

	//---------------------------------------------------------------------------
	// End of Input Section
	//---------------------------------------------------------------------------
		
	// Root stuff
	TFile *file = new TFile(output_filename.c_str(), "RECREATE");
	TTree *VANDMCtree = new TTree("data", "VANDMC output tree");
	
	ReactionProductStructure EJECTdata;
	ReactionProductStructure RECOILdata;
	ReactionObjectStructure REACTIONdata;
	
	VANDMCtree->Branch("eject", &EJECTdata);
	VANDMCtree->Branch("recoil", &RECOILdata);
	if(WriteReaction){
		VANDMCtree->Branch("reaction", &REACTIONdata);
	}

	// Write reaction info to the file.
	std::vector<TNamed*> named;
	SetName(named, "version", VERSION);
	SetName(named, "beamZ", beam_part.GetZ());
	SetName(named, "beamA", beam_part.GetA());
	SetName(named, "beamMass", beam_part.GetMassAMU(), "amu");
	SetName(named, "targZ", targ.GetZ());
	SetName(named, "targA", targ.GetA());
	SetName(named, "targMass", targ.GetMassAMU(), "amu");
	SetName(named, "recoilZ", recoil_part.GetZ());
	SetName(named, "recoilA", recoil_part.GetA());
	SetName(named, "recoilMass", recoil_part.GetMassAMU(), "amu");
	SetName(named, "ejectileZ", eject_part.GetZ());
	SetName(named, "ejectileA", eject_part.GetA());
	SetName(named, "ejectileMass", eject_part.GetMassAMU(), "amu");
	if(InverseKinematics){ SetName(named, "inverse", "Yes"); }
	else{ SetName(named, "inverse", "No"); }
	SetName(named, "beamEnergy", Ebeam0, "MeV");
	SetName(named, "beamType", BeamType);
	SetName(named, "beamspotDiameter", beamspot, "m");
	SetName(named, "beamDivergence", beamAngdiv*rad2deg, "deg");
	SetName(named, "beamEspread", beamEspread, "MeV");
	SetName(named, "Qgs", gsQvalue, "MeV");
	SetName(named, "nExStates", NRecoilStates-1, "MeV");
	SetName(named, "recoilGS", "0.0", "MeV");
	for(unsigned int i = 1; i < NRecoilStates; i++){
		std::stringstream stream; stream << i;
		SetName(named, "recoilState"+stream.str(), ExRecoilStates[i], "MeV");
	}
	if(ADists == 1){ 
		SetName(named, "xsections", "Yes, 1"); 
		for(unsigned int i = 0; i < NRecoilStates; i++){
			std::stringstream stream; stream << i;
			SetName(named, "state"+stream.str()+"Dist", AngDist_fname[i]);
		}
	}
	else if(ADists == 2){
		SetName(named, "xsections", "Yes, 2");
		for(unsigned int i = 0; i < NRecoilStates; i++){
			std::stringstream stream; stream << i;
			SetName(named, "state"+stream.str()+"Rate", AngDist_fname[i], "per event");
		}
	}
	else{ SetName(named, "xsections", "No"); }
	SetName(named, "targetMaterial", targ_mat_name);
	if(targ_mat_name != "NONE"){ SetName(named, "targetMaterial", targ.GetThickness(), "mg/cm^2"); }	
	else{ SetName(named, "targetThickness", targ.GetRealThickness(), "m"); }	
	SetName(named, "targetAngle", targ.GetAngle()*rad2deg, "deg");
	if(PerfectDet){ SetName(named, "perfectDetectors", "Yes"); }
	else{ SetName(named, "perfectDetectors", "No"); }
	SetName(named, "detectorFilename", detector_filename);
	SetName(named, "nDetections", Nwanted);
	if(backgroundRate > 0){ 
		if(bgPerDetection){ SetName(named, "backgroundRate", backgroundRate, "per detection"); }
		else{ SetName(named, "backgroundRate", backgroundRate, "per recoil"); }
		SetName(named, "backgroundWindow", detWindow, "ns");
	}
	else{ SetName(named, "backgroundRate", "NONE"); }
	if(InCoincidence){ SetName(named, "recoilCoincidence", "Yes"); }
	else{ SetName(named, "recoilCoincidence", "No"); }
	if(WriteReaction){ SetName(named, "writeReaction", "Yes"); }
	else{ SetName(named, "writeReaction", "No"); }

	// Create a directory for storing setup information.
	file->mkdir("config");
	file->cd("config");

	// Write the configuration TNameds to file.
	for(std::vector<TNamed*>::iterator iter = named.begin(); iter != named.end(); iter++){
		(*iter)->Write();
		delete (*iter);
	}
	named.clear();
	
	// Create a directory for storing detector setup information.
	file->mkdir("detector");
	file->cd("detector");
	
	// Write all detector entries to file.
	for(size_t index = 0; index < vandle_bars.size(); index++){
		std::stringstream stream; stream << "det";
		if(index < 10){ stream << "0"; }
		stream << index;
		SetName(named, stream.str(), vandle_bars.at(index)->DumpDet());
	}
	
	// Write the detector TNameds to file.
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
	double fpath1 = 0.0, fpath2 = 0.0;
	double recoil_tof = 0.0;
	double eject_tof = 0.0;
	double gamma_tof = 0.0;
	float totTime = 0.0;
	int counter = 1;
	bool flag = false;
	bool hit = false;
	bool veto_event = false;
	unsigned int chunk = Nwanted/10;
	unsigned int beam_stopped = 0;
	unsigned int recoil_stopped = 0;
	unsigned int eject_stopped = 0;
	timer = clock();
	
	int detector_type = 0;
	
	int recoil_detections = 0;
	int eject_detections = 0;
	int gamma_detections = 0;
	
	// Struct for storing reaction information.
	reactData rdata;
	
	while(NgoodDetections < Nwanted){
		// ****************Time Estimate**************
		if(flag && (NgoodDetections % chunk == 0)){
			flag = false;
			totTime = (float)(clock()-timer)/CLOCKS_PER_SEC;
			std::cout << "\n ------------------------------------------------\n"; 
			std::cout << " Number of particles Simulated: " << Nsimulated << std::endl; 
			std::cout << " Number of particles Detected: " << Ndetected << std::endl;
			std::cout << " Number of ejecile particles Detected: " << NgoodDetections << std::endl; 
			if(SupplyRates && ADists){ std::cout << " Number of Reactions: " << Nreactions << std::endl; }
		
			std::cout << " " << NgoodDetections*100.0/Nwanted << "% of simulation complete...\n"; 
			if(PerfectDet){ std::cout << "  Detection Efficiency: " << NgoodDetections*100.0/Nreactions << "%\n"; }
			else{
				std::cout << "  Geometric Efficiency: " << NdetHit*100.0/Nreactions << "%\n";
				std::cout << "  Detection Efficiency: " << NgoodDetections*100.0/Nreactions << "%\n"; 
			}
			if(SupplyRates){ std::cout << "  Beam Time: " << Nsimulated/BeamRate << " seconds\n"; }
		
			std::cout << "  Simulation Time: " << totTime << " seconds\n";
			std::cout << "  Time reamining: " << (totTime/counter)*(10-counter) << " seconds\n";
			counter++; 
		}

		recoil_detections = 0;
		eject_detections = 0;
		gamma_detections = 0;

		if(backgroundWait != 0){ // Simulating background events
			backgroundWait--;

			// Process the background event for each detector
			for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
				if(!((*iter)->IsEjectileDet() || (*iter)->IsRecoilDet())){ continue; } // This detector cannot detect particles

				// Select the "tof" of the background event. Use the recoil tof, because it doesn't matter.
				recoil_tof = (double)frand(0, detWindow)*(1E-9);
			
				// Select a random point insde the detector
				(*iter)->GetRandomPointInside(temp_vector);
				Cart2Sphere(temp_vector_sphere);
			
				// Calculate the apparent energy of the particle using the tof
				if((*iter)->IsEjectileDet()){
					EJECTdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], temp_vector.Length(), temp_vector_sphere.axis[1]*rad2deg,
									 temp_vector_sphere.axis[2]*rad2deg, 0.0, recoil_tof*(1E9), 0.0, 0.0, 0.0, 0.0, (*iter)->GetLoc(), true);
					VANDMCtree->Fill(); 
					EJECTdata.Zero();
				}
				else if((*iter)->IsRecoilDet()){
					RECOILdata.Append(temp_vector.axis[0], temp_vector.axis[1], temp_vector.axis[2], temp_vector.Length(), RecoilSphere.axis[1]*rad2deg,
									  RecoilSphere.axis[2]*rad2deg, 0.0, recoil_tof*(1E9), 0.0, 0.0, 0.0, 0.0, (*iter)->GetLoc(), true);
					VANDMCtree->Fill();
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
			if(lab_beam_focus.axis[2] != 0.0){ 
				// In this case, lab_beam_focus is the originating point of the beam particle
				// (or the terminating point for beams focused downstream of the target)
				// The direction is given by the cartesian vector 'lab_beam_trajectory'
				if(BeamType == 0){ RandomCircle(beamspot, lab_beam_focus.axis[2], lab_beam_trajectory); } // Gaussian beam
				else if(BeamType == 1){ RandomGauss(beamspot, lab_beam_focus.axis[2], lab_beam_trajectory); } // Cylindrical beam
				else if(BeamType == 2){ RandomHalo(beamspot/2.0, lab_beam_focus.axis[2], lab_beam_trajectory); } // Halo beam
				
				// Normalize the trajectory and calculate the interaction depth.
				lab_beam_trajectory.Normalize();
				Zdepth = targ.GetInteractionDepth(lab_beam_focus, lab_beam_trajectory, targ_surface, lab_beam_interaction);
			}
			else{ 
				// In this case, lab_beam_start stores the originating point of the beam particle
				// The direction is given simply by the +z-axis
				// The 1m offset ensures the particle originates outside the target
				if(BeamType == 0){ RandomGauss(beamspot/2.0, -1.0, lab_beam_start); } // Gaussian beam
				else if(BeamType == 1){ RandomCircle(beamspot, -1.0, lab_beam_start); } // Cylindrical beam
				else if(BeamType == 2){ RandomHalo(beamspot/2.0, -1.0, lab_beam_start); } // Halo beam
				lab_beam_start.axis[2] *= -1; // This is done to place the beam particle upstream of the target
				
				// Normalize the trajectory and calculate the interaction depth.
				lab_beam_trajectory.Normalize();
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
			Egamma = rdata.Eexcited;

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
		}

		// Set the first detector type to process.
		if(have_veto_det){ detector_type = -1; }
		else if(have_recoil_det){ detector_type = 0; }
		else if(have_ejectile_det){ detector_type = 1; }
		else if(have_gamma_det){ detector_type = 2; }

		recoil_tof = -1;
		veto_event = false;

process:
		// Process the reaction products
		for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
			// Check if we need to process this detector.
			if(detector_type == -1){
				if(!(*iter)->IsVeto()){ continue; } // This detector is not a veto detector.
			}
			else if(detector_type == 0){
				if(!(*iter)->IsRecoilDet()){ continue; } // This detector is not able to detect recoils.
				if(ErecoilMod <= 0.0){ break; } // The recoil has stopped. We're done tracking it.
			}
			else if(detector_type == 1){
				if(!(*iter)->IsEjectileDet()){ continue; } // This detector is not able to detect ejectiles.
				if(EejectMod <= 0.0){ break; } // The ejectile has stopped. We're done tracking it.
			}
			else if(detector_type == 2){
				if(!(*iter)->IsGammaDet()){ continue; } // This detector is not able to detect gammas.
				if(Egamma <= 0.0 && !NeutronSource){ break; } // Do not process the ground state.
			}
			else{ continue; } // This detector cannot detect particles.

			if(detector_type == -1){ // Process the ejectile and recoil for veto detectors.
				veto_event = (*iter)->IntersectPrimitive(lab_beam_interaction, Recoil, HitDetect1, fpath1, fpath2) ||
				             (*iter)->IntersectPrimitive(lab_beam_interaction, Ejectile, HitDetect1, fpath1, fpath2); 
				      
				// If the veto detector was triggered, stop this event.
				if(veto_event){ break; }
			}
			else if(detector_type == 0){ // Process the recoil.
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
			else if(detector_type == 1){ // Process the ejectile.
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
			else if(detector_type == 2){ // Process the gamma ray.
				// Simulate the gamma emission. The gamma rays are emitted isotropically from the reaction point.
				UnitSphereRandom(dummy_vector);
				hit = (*iter)->IntersectPrimitive(lab_beam_interaction, dummy_vector, HitDetect1, fpath1, fpath2);
			}
			
			// If a geometric hit was detected, process the particle
			if(hit){
				NdetHit++; 

				// Solve for the energy deposited in the material.
				if((*iter)->UseMaterial()){ // Do energy loss and range considerations
					if(detector_type == 0){ 
						if(recoil_part.GetZ() > 0){ // Calculate energy loss for the recoil in the detector
							QDC = ErecoilMod - recoil_tables[(*iter)->GetMaterial()].GetNewE(ErecoilMod, temp_vector.Length(), dist_traveled);
						}
						else{ std::cout << " ERROR: Doing energy loss on recoil particle with Z == 0???\n"; }
					}	
					else if(detector_type == 1){
						if(eject_part.GetZ() > 0){ // Calculate energy loss for the ejectile in the detector
							QDC = EejectMod - eject_tables[(*iter)->GetMaterial()].GetNewE(EejectMod, temp_vector.Length(), dist_traveled);
						}
						else{ std::cout << " ERROR: Doing energy loss on ejectile particle with Z == 0???\n"; }
					}
					else if(detector_type == 2){ std::cout << " ERROR: Doing energy loss on a gamma ray???\n"; }
				}
				else{ // Do not do energy loss calculations. The particle leaves all of its energy in the detector.
					dist_traveled = temp_vector.Length()*frand(); // The particle penetrates a random distance into the detector and stops.
					if(detector_type == 0){ QDC = ErecoilMod; } // The recoil may leave any portion of its energy inside the detector
					else if(detector_type == 1){ QDC = EejectMod; } // The ejectile may leave any portion of its energy inside the detector
					else if(detector_type == 2){ QDC = Egamma; }
				}

				// If particle originates outside of the detector, add the flight path to the first encountered
				// detector face. Otherwise, if the particle originates inside the detector (i.e. a detector at
				// the origin), the distance traveled through the detector is already the total flight path and
				// the time of flight is irrelevant.
				if(fpath1 >= 0.0){
					dist_traveled += fpath1;
					
					// Calculate the particle ToF (ns)
					if(detector_type == 0){ 
						recoil_tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMrecoilMeV()/ErecoilMod); 
						recoil_tof += rndgauss0(timeRes); // Smear tof due to PIXIE resolution
					}
					else if(detector_type == 1){ 
						eject_tof = (dist_traveled/c)*std::sqrt(0.5*kind.GetMejectMeV()/EejectMod);
						eject_tof += rndgauss0(timeRes); // Smear tof due to PIXIE resolution
						if(recoil_tof > 0.0){ eject_tof = eject_tof - recoil_tof; }
					}
					else if(detector_type == 2){ 
						gamma_tof = dist_traveled/c;
						gamma_tof += rndgauss0(timeRes); // Smear tof due to PIXIE resolution
						if(recoil_tof > 0.0){ gamma_tof = gamma_tof - recoil_tof; }
					}
				}

				// Get the local coordinates of the intersection point.
				(*iter)->GetLocalCoords(HitDetect1, hit_x, hit_y, hit_z);
			
				// Calculate the lab angles of the detector intersection point. Ignore normalization, we're
				// going to throw away R anyway.
				if(detector_type == 0){ Cart2Sphere(HitDetect1, RecoilSphere); }
				else if(detector_type == 1){ Cart2Sphere(HitDetect1, EjectSphere); }
				else if(detector_type == 2){ Cart2Sphere(HitDetect1, GammaSphere); }
			
				// Calculate the hit detection point in 3d space. This point will lie along the vector pointing
				// from the origin to the point where the ray intersects a detector and takes finite range of a
				// particle in a material into account.
				HitDetect1.Normalize();
				HitDetect1 = HitDetect1*dist_traveled;
			
				// Main output
				if(detector_type == 0){
					RECOILdata.Append(HitDetect1.axis[0], HitDetect1.axis[1], HitDetect1.axis[2], HitDetect1.Length(), RecoilSphere.axis[1]*rad2deg,
					                  RecoilSphere.axis[2]*rad2deg, QDC, recoil_tof*(1E9), rdata.Erecoil, hit_x, hit_y, hit_z, (*iter)->GetLoc(), false);
				
					recoil_detections++;
				
					// Adjust the recoil energy to take energy loss into account. 
					ErecoilMod = ErecoilMod - QDC;
				}
				else if(detector_type == 1){
					EJECTdata.Append(HitDetect1.axis[0], HitDetect1.axis[1], HitDetect1.axis[2], HitDetect1.Length(), EjectSphere.axis[1]*rad2deg,
					                  EjectSphere.axis[2]*rad2deg, QDC, eject_tof*(1E9), rdata.Eeject, hit_x, hit_y, hit_z, (*iter)->GetLoc(), false);
								
					eject_detections++;			
								
					// Adjust the ejectile energy to take energy loss into account. 
					EejectMod = EejectMod - QDC;
				}
				else if(detector_type == 2){ 
					EJECTdata.Append(HitDetect1.axis[0], HitDetect1.axis[1], HitDetect1.axis[2], HitDetect1.Length(), GammaSphere.axis[1]*rad2deg,
					                 GammaSphere.axis[2]*rad2deg, Egamma, gamma_tof*(1E9), 0.0, hit_x, hit_y, hit_z, (*iter)->GetLoc(), true);
									 
					gamma_detections++;
					
					// Done tracking the gamma ray.		 
					break;
				}
			} // if(hit)
		} // for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++)
		
		if(!veto_event){ 
			// Decide which detector type to process next, if any.
			if(detector_type == -1){
				if(have_recoil_det){
					detector_type = 0;
					goto process;
				}
				else if(have_ejectile_det){ 
					detector_type = 1; 
					goto process;
				}
				else if(have_gamma_det){
					detector_type = 2;
					goto process;
				}
			}
			else if(detector_type == 0){
				if(have_ejectile_det){ 
					detector_type = 1; 
					goto process;
				}
				else if(have_gamma_det){
					detector_type = 2;
					goto process;
				}
			}
			else if(detector_type == 1){
				if(have_gamma_det){
					detector_type = 2;
					goto process;
				}
			}
			else if(detector_type == 2){
			}

			NrecoilHits += recoil_detections;
			NejectileHits += eject_detections;
			NgammaHits += gamma_detections;

			// Check to see if anything needs to be written to file.
			if(InCoincidence){ // We require coincidence between ejectiles and recoils 
				if(recoil_detections > 0 && (eject_detections > 0 || gamma_detections > 0)){ 
					if(WriteReaction){ // Set some extra reaction data variables.
						REACTIONdata.Append(rdata.Ereact, rdata.Eeject, rdata.Erecoil, rdata.comAngle*rad2deg, rdata.state,
							                lab_beam_interaction.axis[0], lab_beam_interaction.axis[1], lab_beam_interaction.axis[2],
							                lab_beam_stragtraject.axis[0], lab_beam_stragtraject.axis[1], lab_beam_stragtraject.axis[2]);
					}
					if(bgPerDetection){ backgroundWait = backgroundRate; }
					VANDMCtree->Fill(); 
					Ndetected++;
					
					// Ignore background and gamma events for true detection count.
					if(eject_detections > 0){ 
						if(!flag){ flag = true; }
						NgoodDetections++; 
					}
				}
			}
			else{ // Coincidence is not required between reaction particles
				if(eject_detections > 0 || recoil_detections > 0 || gamma_detections > 0){ 
					if(WriteReaction){
						REACTIONdata.Append(rdata.Ereact, rdata.Eeject, rdata.Erecoil, rdata.comAngle*rad2deg, rdata.state,
							                lab_beam_interaction.axis[0], lab_beam_interaction.axis[1], lab_beam_interaction.axis[2],
							                lab_beam_stragtraject.axis[0], lab_beam_stragtraject.axis[1], lab_beam_stragtraject.axis[2]);
					}
					if(bgPerDetection){ backgroundWait = backgroundRate; }
					VANDMCtree->Fill(); 
					Ndetected++;
					
					// Ignore background and gamma events for true detection count.
					if(eject_detections > 0){ 
						if(!flag){ flag = true; }
						NgoodDetections++; 
					}
				}
			}
		}
		else{ NvetoEvents++; }
		
		// Zero all output data structures.
		EJECTdata.Zero();
		RECOILdata.Zero();
		if(WriteReaction){ REACTIONdata.Zero(); }
	} // Main simulation loop
	// ==  ==  ==  ==  ==  ==  == 

	// Create a directory for storing end of simulation information.
	file->mkdir("simulation");
	file->cd("simulation");

	SetName(named, "simulationTime", (float)(clock()-timer)/CLOCKS_PER_SEC, "seconds");
	SetName(named, "totalEvents", Nreactions);
	SetName(named, "totalDetectorHits", NdetHit);
	SetName(named, "totalDetectedEvents", Ndetected);
	SetName(named, "vetoedEvents", NvetoEvents);
	SetName(named, "recoilHits", NrecoilHits);
	SetName(named, "ejectileHits", NejectileHits);
	SetName(named, "gammaHits", NgammaHits);

	// Write the configuration TNameds to file.
	for(std::vector<TNamed*>::iterator iter = named.begin(); iter != named.end(); iter++){
		(*iter)->Write();
		delete (*iter);
	}
	named.clear();
	
	// Information output and cleanup
	std::cout << "\n ------------- Simulation Complete --------------\n";
	std::cout << " Simulation Time: " << (float)(clock()-timer)/CLOCKS_PER_SEC << " seconds\n"; 
	std::cout << " Total MC Events: " << Nreactions << "\n";
	std::cout << " Total Detector Hits: " << NdetHit << "\n";
	std::cout << "  Vetoed Events: " << NvetoEvents << " (" << (100.0*NvetoEvents)/Nreactions << "%)\n";
	std::cout << "  Recoil Hits:   " << NrecoilHits << " (" << (100.0*NrecoilHits)/Nreactions << "%)\n";
	std::cout << "  Ejectile Hits: " << NejectileHits << " (" << (100.0*NejectileHits)/Nreactions << "%)\n";
	std::cout << "  Gamma Hits:    " << NgammaHits << " (" << (100.0*NgammaHits)/Nreactions << "%)\n";
	if(beam_stopped > 0 || eject_stopped > 0 || recoil_stopped > 0){
		std::cout << " Particles Stopped in Target:\n";
		if(beam_stopped > 0){ std::cout << "  Beam: " << beam_stopped << " (" << 100.0*beam_stopped/Nsimulated << "%)\n"; }
		if(eject_stopped > 0){ std::cout << "  Ejectiles: " << eject_stopped << " (" << 100.0*eject_stopped/Nsimulated << "%)\n"; }
		if(recoil_stopped > 0){ std::cout << "  Recoils: " << recoil_stopped << " (" << 100.0*recoil_stopped/Nsimulated << "%)\n"; }
	}
	if(SupplyRates){ 
		double beamTime = 0.0;
		for(unsigned int i = 0; i < NRecoilStates; i++){
			beamTime += kind.GetDistribution(i)->GetRate()*kind.GetNreactions(i);
		}
		std::cout << " Beam Time: " << beamTime << " seconds (" << beamTime/3600 << " hrs.)\n"; 
	}
	
	file->cd();
	VANDMCtree->Write();

	std::cout << "  Wrote file " << output_filename << "\n";
	std::cout << "   Wrote " << VANDMCtree->GetEntries() << " tree entries for VANDMC\n";
	file->Close();
	
	delete file;
	delete[] ExRecoilStates;
	delete[] totXsect;

	for(std::vector<Primitive*>::iterator iter = vandle_bars.begin(); iter != vandle_bars.end(); iter++){
		delete *iter;
	}
	vandle_bars.clear();
	
	return true;
} 

int main(int argc, char *argv[]){
	vandmc obj;
	return (obj.Execute(argc, argv) ? 0 : 1);
}
