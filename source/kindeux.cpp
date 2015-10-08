// kindeux.cpp
// Cory Thornsberry

#include "vikar_core.h"
#include "kindeux.h"

/////////////////////////////////////////////////////////////////////
// Kindeux
/////////////////////////////////////////////////////////////////////
    	
// Mass values are input as AMU
// The ground state Q-value is given in MeV
// NrecoilStates_ is the number of excited states of the recoil 
// RecoilExStates_ is a pointer to an array of excitations for the recoil (in MeV)
void Kindeux::Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, unsigned int NrecoilStates_, double *RecoilExStates_){
	if(!init){
		Mbeam = Mbeam_; 
		Mtarg = Mtarg_; 
		Mrecoil = Mrecoil_; 
		Meject = Meject_; 
		Qvalue = Qvalue_;
		NrecoilStates = NrecoilStates_;
		RecoilExStates = RecoilExStates_;
		init = true;
	}
}

// Set Kindeux to use angular distributions for calculating ejectile angles
// Returns false if attempt to load the distributions fails for any reason
// tgt_thickness_ is given in units of mg/cm^2
bool Kindeux::SetDist(std::vector<std::string> &fnames, double total_targ_mass, double tgt_thickness_, double incident_beam_current){
	if(!init || ang_dist){ return false; }
	if(fnames.size() < NrecoilStates){
		std::cout << " Kindeux: Warning! Must have distributions for " << NrecoilStates << " excited states and the ground state\n";
		std::cout << " Kindeux: Received distributions for only " << fnames.size() << " states but expected " << NrecoilStates << std::endl;
		return false;
	}
	
	ang_dist = true;
	distributions = new AngularDist[NrecoilStates];
	Xsections = new double[NrecoilStates];

	// Load all distributions from file
	total_xsection = 0.0;
	for(unsigned int i = 0; i < NrecoilStates; i++){
		if(!distributions[i].Initialize(fnames[i].c_str(), total_targ_mass, tgt_thickness_, incident_beam_current)){
			std::cout << "  Failed to load angular distribution file '" << fnames[i] << "'\n";
			ang_dist = false;
			break;
		}
		Xsections[i] = total_xsection;
		total_xsection += distributions[i].GetReactionXsection();
	}

	// Encountered some problem with one or more of the distributions
	if(!ang_dist){ 
		delete[] distributions; 
		delete[] Xsections;
		return false;
	}
	return true;
}

// Set Kindeux to use Rutherford scattering for calculating ejectile angles
// coeff_ in m^2
bool Kindeux::SetRutherford(double coeff_){
	if(!init || ang_dist){ return false; }
	
	ang_dist = true;
	distributions = new AngularDist[1];
	Xsections = new double[1]; Xsections[0] = 0.0;

	double angles[181], xsections[181];
	for(unsigned int i = 1; i <= 180; i++){
		angles[i] = (double)i; // degrees
		xsections[i] = coeff_*(1.0 / pow(std::sin((angles[i]/2.0)*deg2rad), 4.0))*(1E28); // mb/sr
	}
	angles[0] = 0; 
	xsections[0] = xsections[1];
	
	distributions[0].Initialize(181, angles, xsections);
	total_xsection = distributions[0].GetReactionXsection();
	
	return true;
}

// Calculate recoil excitation energies
// Returns true if there was a reaction and false otherwise
bool Kindeux::get_excitations(double &recoilE, unsigned int &state){
	if(NrecoilStates == 0){
		state = 0;
		recoilE = RecoilExStates[state];
		return true;
	}
	else if(ang_dist){	
		// Angular dist weighted
		double rand_xsection = frand()*total_xsection;
		for(unsigned int i = 0; i < NrecoilStates-1; i++){
			if(rand_xsection >= Xsections[i] && rand_xsection <= Xsections[i+1]){
				// State i has been selected
				state = i;
				recoilE = RecoilExStates[state];
				return true;
			}
		}
		// rand_xsection falls in the range (Xsections[NrecoilStates-1], total_xsection]
		// State NrecoilStates-1 has been selected
		state = NrecoilStates-1;
		recoilE = RecoilExStates[state];
		return true;
	}
	else{ 
		// Isotropic
		state = (unsigned int)(frand()*NrecoilStates);
		recoilE = RecoilExStates[state];
	}
	return true;
}

// See J. B. Ball, "Kinematics II: A Non-Relativistic Kinematics Fortran Program
// to Aid Analysis of Nuclear Reaction Angular Distribution Data", ORNL-3251
bool Kindeux::FillVars(reactData &react, Vector3 &Ejectile, int recoil_state/*=-1*/, int solution/*=-1*/, double theta/*=-1*/){
	Vector3 dummy_Recoil;
	return FillVars(react, Ejectile, dummy_Recoil, recoil_state, solution, theta);
}

// Overloaded version which also calculates data for the recoil particle
bool Kindeux::FillVars(reactData &react, Vector3 &Ejectile, Vector3 &Recoil, int recoil_state/*=-1*/, int solution/*=-1*/, double theta/*=-1*/){
	double Recoil_Ex;
	if(recoil_state >= 0 && (unsigned int)recoil_state < NrecoilStates){ // Select the state to use
		react.state = recoil_state;
		Recoil_Ex = RecoilExStates[react.state];
	}
	else if(!get_excitations(Recoil_Ex, react.state)){ // Use built-in state selection
		// No reaction occured
		return false; 
	}
	
	// In the center of mass frame
	double EjectPhi, EjectTheta;
	double Vcm = std::sqrt(2.0*Mbeam*react.Ereact)/(Mbeam+Mtarg); // Velocity of the center of mass
	double Ecm = Mtarg*react.Ereact/(Mbeam+Mtarg); // Energy of the center of mass

	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	react.comAngle = -1.0; // Ejectile and recoil angle in the center of mass frame
	
	if(theta >= 0.0){
		double temp;
		UnitSphereRandom(temp, EjectPhi);
		
		if(theta >= 0.0){ react.comAngle = theta; }
		else{ react.comAngle = temp; }
	}
	else if(ang_dist){
		// Sample the angular distributions for the CoM angle of the ejectile
		if(distributions[react.state].Sample(react.comAngle)){ EjectPhi = 2*pi*frand(); } // Randomly select phi of the ejectile
		else{ UnitSphereRandom(react.comAngle, EjectPhi); } // Failed to sample the distribution
	}
	else{ UnitSphereRandom(react.comAngle, EjectPhi); } // Randomly select a uniformly distributed point on the unit sphere
	
	EjectTheta = std::atan2(std::sin(react.comAngle),(std::cos(react.comAngle)+(Vcm/VejectCoM))); // Ejectile angle in the lab
	double temp_value = std::sqrt(VejectCoM*VejectCoM-pow(Vcm*std::sin(EjectTheta),2.0));
	double Ejectile_V = Vcm*std::cos(EjectTheta); // Ejectile velocity in the lab frame
	
	if(VejectCoM >= Vcm){ 
		// Veject is single valued
		Ejectile_V += temp_value; 
	} 
	else{ 
		// Veject is double valued, so we randomly choose one of the values
		// for the velocity, and hence, the energy of the ejectile
		if(solution < 0){
			if(frand() >= 0.5){ Ejectile_V += temp_value; }
			else{ Ejectile_V = Ejectile_V - temp_value; }
		}
		else if(solution == 0){ Ejectile_V += temp_value; }
		else{ Ejectile_V = Ejectile_V - temp_value; }
	}
	
	react.Eeject = 0.5*Meject*Ejectile_V*Ejectile_V; // Ejectile energy in the lab frame
	Ejectile = Vector3(1.0, EjectTheta, EjectPhi); // Ejectile direction unit vector
	
	// Recoil calculations (now in the lab frame)
	react.Erecoil = (react.Ereact+Qvalue-(0.0+Recoil_Ex)) - react.Eeject;
	Recoil = Vector3(1.0, std::asin(((std::sqrt(2*Meject*react.Eeject))/(std::sqrt(2*Mrecoil*react.Erecoil)))*std::sin(EjectTheta)), WrapValue(EjectPhi+pi,0.0,2*pi));
	
	return true;
}

// Convert ejectile CoM angle to Lab angle
double Kindeux::ConvertAngle2Lab(double Beam_E, double Recoil_Ex, double Eject_CoM_angle){
	// In the center of mass frame
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass	
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	
	return(std::atan2(std::sin(Eject_CoM_angle),(std::cos(Eject_CoM_angle)+(Vcm/VejectCoM)))); // Ejectile angle in the lab
}

// Print debug information
// Does nothing if angular distributions are not set
void Kindeux::Print(){
	if(ang_dist){
		for(unsigned int i = 0; i < NrecoilStates; i++){
			std::cout << "  State " << i+1 << ":";
			std::cout << " Reaction X-Section: " << distributions[i].GetReactionXsection() << " mb";
			std::cout << "\tExpected Rate: " << distributions[i].GetRate() << " pps\n";
		}
	}
}
