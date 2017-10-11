/** \file kindeux.cpp
 * \brief Performs kinematics calculations.
 *
 * The Kindeux class is intended to be used for calculating energies
 * and outgoing angles for two body reactions. The kinematics methods
 * are based on a document by J. B. Ball (ORNL-3251) and are used in
 * conjunction with a monte carlo program to simulate low energy
 * reactions.
 *
 * \author C. R. Thornsberry
 * \date Feb. 26th, 2016
 */
#ifndef KINDEUX_H
#define KINDEUX_H

#include <vector>
#include <string>

class AngularDist;
class Target;

class reactData{
  public:
	double Ereact;
	double Eeject;
	double Erecoil;
	double Eexcited;
	double comAngle;
	unsigned int state;

	reactData() : Ereact(0.0), Eeject(0.0), Erecoil(0.0), Eexcited(0.0), comAngle(0.0), state(0) { }
};

class Californium{
  private:
        double energy[101];
        double neutrons[101];
        double integral[101];
        double totalIntegral;

        const double SFbranch = 0.0309; 
        const double AlphaBranch = 0.9691;

        double func(const double &E);

  public:
        Californium();

        double sample();
};

class Kindeux{
  private:   	
	double Mbeam, Mtarg, Mrecoil;
	double Meject, Qvalue;
	double *RecoilExStates;
	double *Xsections;
	double total_xsection;
	
	unsigned int *Nreactions;
	unsigned int NDist, NrecoilStates;
	AngularDist *distributions;
	bool ang_dist, init;
	bool inverse;
	bool nsource;
	
	Californium cf;

	/// Get the excitation of the recoil particle.
	bool get_excitations(double &recoilE, unsigned int &state);
	
  public:
	Kindeux();
	
	~Kindeux();

	/// Return true if the object is initialized, and false otherwise.
	bool IsInit(){ return init; }

	/// Return the mass of the beam particle in amu or kg.
	double GetMbeam(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } 
	
	/// Return the mass of the target particle in amu or kg.
	double GetMtarg(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } 
	
	/// Return the mass of the recoil particle in amu or kg.
	double GetMrecoil(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } 
	
	/// Return the mass of the ejectile particle in amu or kg.
	double GetMeject(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } 

	/// Return the mass of the beam particle in MeV/c^2.
	double GetMbeamMeV(){ return 931.49*Mbeam; } 

	/// Return the mass of the target particle in MeV/c^2.
	double GetMtargMeV(){ return 931.49*Mtarg; } 

	/// Return the mass of the recoil particle in MeV/c^2.
	double GetMrecoilMeV(){ return 931.49*Mrecoil; } 

	/// Return the mass of the ejectile particle in MeV/c^2.
	double GetMejectMeV(){ return 931.49*Meject; } 
   	
   	/// Return the number of recoil states being used.
   	unsigned int GetNrecoilStates(){ return NrecoilStates; }
   	
   	/// Return the number of reactions for a given state.
   	unsigned int GetNreactions(const unsigned int &index_){ return (index_ < NrecoilStates ? Nreactions[index_] : 0); }
   	
   	/// Return a pointer to the angular distribution for a given state.
   	AngularDist *GetDistribution(const unsigned int &index_){ return (ang_dist && index_ < NrecoilStates ? &distributions[index_] : NULL); }
   
	/// Return true if this reaction is in inverse kinematics and false otherwise.
	bool GetInverseKinematics(){ return inverse; }
	
   	/// Initialize Kindeux object with reaction parameters.
   	void Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, unsigned int NrecoilStates_, double *RecoilExStates_);
	
	/// Set Kindeux to use angular distributions from files for calculating recoil excitations.
	bool SetDist(std::vector<std::string> &fnames, double incident_beam_current, Target *targ_);
	
	/// Set Kindeux to use relative state intensities for calculating recoil excitations.
	bool SetDist(const std::vector<std::string> &intensities_);
	
	/// Set Kindeux to use Rutherford scattering for calculating reaction product angles.
	bool SetRutherford(double coeff_);
	
	/// Calculate reaction product energy and angle for the ejectile particle only.
	bool FillVars(reactData &react, Vector3 &Ejectile, int recoil_state=-1, int solution=-1, double theta=-1);
	
	/// Calculate reaction product energies and angles for the recoil and ejectile particles.
	bool FillVars(reactData &react, Vector3 &Ejectile, Vector3 &Recoil, int recoil_state=-1, int solution=-1, double theta=-1);
	
	/// Toggle whether or not to use this class as a neutron source.
	bool ToggleNeutronSource(){ return (nsource = !nsource); }

	/// Convert an input center of mass angle to the lab frame.
	double ConvertAngle2Lab(double, double, double);
	
	/// Print information about the kindeux reaction object.
	void Print();
};

#endif
