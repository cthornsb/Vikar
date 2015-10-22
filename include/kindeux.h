// kindeux.h
// Cory Thornsberry

#ifndef KINDEUX_H
#define KINDEUX_H

#include <vector>
#include <string>

class AngularDist;

struct reactData{
	double Ereact;
	double Eeject;
	double Erecoil;
	double comAngle;
	unsigned int state;
};

class Kindeux{
  private:   	
	double Mbeam, Mtarg, Mrecoil;
	double Meject, Qvalue;
	double *RecoilExStates;
	double *Xsections;
	double total_xsection;
	
	unsigned int NDist, NrecoilStates;
	AngularDist *distributions;
	bool ang_dist, init;
	
	bool get_excitations(double &recoilE, unsigned int &state);
	
  public:
	Kindeux(){
		ang_dist = false; init = false;
		NDist = 0; NrecoilStates = 0;
		RecoilExStates = NULL;
		Mbeam = 0.0; Mtarg = 0.0;
		Mrecoil = 0.0; Meject = 0.0;
		Qvalue = 0.0;
	}
	~Kindeux(){ 
		if(ang_dist){ 
			delete[] distributions; 
			delete[] Xsections;
		} 
	}

	bool IsInit(){ return init; }

	double GetMbeam(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } /// Return the mass of the beam particle in amu or kg
	double GetMtarg(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } /// Return the mass of the beam particle in amu or kg
	double GetMrecoil(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } /// Return the mass of the beam particle in amu or kg
	double GetMeject(bool in_kg=false){ return (!in_kg ? Mbeam : Mbeam/6.02214129E26); } /// Return the mass of the beam particle in amu or kg

	double GetMbeamMeV(){ return 931.49*Mbeam; } /// Return the mass of the beam particle in MeV/c^2
	double GetMtargMeV(){ return 931.49*Mtarg; } /// Return the mass of the target particle in MeV/c^2
	double GetMrecoilMeV(){ return 931.49*Mrecoil; } /// Return the mass of the recoil particle in MeV/c^2
	double GetMejectMeV(){ return 931.49*Meject; } /// Return the mass of the ejectile particle in MeV/c^2
   	
   	void Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, unsigned int NrecoilStates_, double *RecoilExStates_);
	bool SetDist(std::vector<std::string> &fnames, double total_targ_mass, double tgt_thickness_, double incident_beam_current);
	bool SetRutherford(double coeff_);
	bool FillVars(reactData &react, Vector3 &Ejectile, int recoil_state=-1, int solution=-1, double theta=-1);
	bool FillVars(reactData &react, Vector3 &Ejectile, Vector3 &Recoil, int recoil_state=-1, int solution=-1, double theta=-1);
	double ConvertAngle2Lab(double, double, double);
	void Print();
};

#endif