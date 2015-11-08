// kindeux.cpp
// Cory Thornsberry

#include "vikar_core.h"
#include "kindeux.h"

/////////////////////////////////////////////////////////////////////
// Kindeux
/////////////////////////////////////////////////////////////////////
    	
/** Initialize Kindeux object with reaction parameters.
  * param[in] Mbeam_ Mass of the beam (amu).
  * param[in] Mtarg_ Mass of the target (amu).
  * param[in] Mrecoil_ Mass of the recoil particle (amu).
  * param[in] Meject_ Mass of the ejectile (amu).
  * param[in] Qvalue_ The ground state Q-value for the reaciton (MeV).
  * param[in] NrecoilStates_ The number of recoil particle excitations.
  * param[in] RecoilExStates_ The array containing all recoil excitations (MeV).
  */
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

/** Set Kindeux to use angular distributions from files for calculating recoil excitations.
  * Returns false if attempt to load distribution fails for any reason.
  * param[in] fnames_ Vector containing filenames for each recoil state distribution (including g.s.).
  * param[in] total_targ_mass_ Molar mass of the target (g/mol).
  * param[in] tgt_thickness_ The thickness of the target (mg/cm^2).
  * param[in] incident_beam_current The incident rate of the beam (pps).
  */
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

/** Set Kindeux to use relative state intensities for calculating recoil excitations.
  * Returns false if attempt to set distribution fails for any reason.
  * param[in] intensities_ Vector containing relative intensities for each recoil state (including g.s.).
  */
bool Kindeux::SetDist(const std::vector<std::string> &intensities_){
	if(!init || ang_dist){ return false; }
	
	if(intensities_.size() < NrecoilStates){
		std::cout << " Kindeux: Warning! Must have state intensities for " << NrecoilStates << " excited states and the ground state\n";
		std::cout << " Kindeux: Received intensities for only " << intensities_.size() << " states but expected " << NrecoilStates << std::endl;
		return false;
	}
	
	ang_dist = true;
	
	distributions = new AngularDist[NrecoilStates];
	Xsections = new double[NrecoilStates];

	// Load all distributions from file
	total_xsection = 0.0;
	for(unsigned int i = 0; i < NrecoilStates; i++){
		if(!distributions[i].Initialize((double)atof(intensities_[i].c_str()))){
			std::cout << "  Failed to set relative intensity for state " << i << std::endl;
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

/** Set Kindeux to use Rutherford scattering for calculating reaction product angles.
  * Returns false if the object has already been initialized.
  * param[in] coeff_ The Rutherford coefficient to use for the distribution (m^2).
  */
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

/** Get the excitation of the recoil particle based on angular distributions
  *  if they are available or isotropic distributions if they are not.
  * Returns true if a reaction occured and false otherwise.
  * param[out] recoilE The excitation energy of the recoil particle.
  * param[out] state The excitation state of the recoil particle.
  */
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

/** Calculate reaction product energy and angle for the ejectile particle only.
  */
bool Kindeux::FillVars(reactData &react, Vector3 &Ejectile, int recoil_state/*=-1*/, int solution/*=-1*/, double theta/*=-1*/){
	Vector3 dummy_Recoil;
	return FillVars(react, Ejectile, dummy_Recoil, recoil_state, solution, theta);
}

/** Calculate reaction product energies and angles for the recoil and ejectile particles.
  * Returns false if no reaction event occured due to the angular distribution and beam rate.
  * See J. B. Ball, "Kinematics II: A Non-Relativistic Kinematics Fortran Program
  *  to Aid Analysis of Nuclear Reaction Angular Distribution Data", ORNL-3251.
  * param[out] react Reaction data object for storing reaction event information.
  * param[out] Ejectile The spherical coordinates vector of the ejectile particle in the lab frame.
  * param[out] Recoil The spherical coordinates vector of the recoil particle in the lab frame.
  * param[in] recoil_state Allows the user to specify the recoil excited state which to use for the reaction.
  * param[in] solution Allows the user to select either the + or - solution of the ejectile velocity.
  * param[in] theta Allows the user to select the reaction center of mass angle (in rad).
  */
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
		react.comAngle = distributions[react.state].Sample();
		if(react.comAngle > 0.0){ EjectPhi = 2*pi*frand(); } // Randomly select phi of the ejectile
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

/** Convert an input center of mass angle to the lab frame.
  * Return the lab frame angle (in rad).
  * param[in] Beam_E the energy of the beam particle (MeV).
  * param[in] Recoil_Ex the excitation of the outgoing recoil particle (MeV).
  * param[in] Eject_CoM_angle the reaction center of mass angle (rad).
  */
double Kindeux::ConvertAngle2Lab(double Beam_E, double Recoil_Ex, double Eject_CoM_angle){
	// In the center of mass frame
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass	
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	
	return(std::atan2(std::sin(Eject_CoM_angle),(std::cos(Eject_CoM_angle)+(Vcm/VejectCoM)))); // Ejectile angle in the lab
}

/// Print information about the kindeux reaction object.
void Kindeux::Print(){
	if(ang_dist){
		for(unsigned int i = 0; i < NrecoilStates; i++){
			std::cout << "  State " << i+1 << ":";
			std::cout << " Reaction X-Section: " << distributions[i].GetReactionXsection() << " mb";
			std::cout << "\tExpected Rate: " << distributions[i].GetRate() << " pps\n";
		}
	}
}
