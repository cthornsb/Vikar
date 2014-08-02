// vikar_lib.h
// Cory Thornsberry

#ifndef VIKAR_LIB_H
#define VIKAR_LIB_H

#include <cmath>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

extern const double c, pi, deg2rad;
extern const double rad2deg, LN2, avagadro;

/////////////////////////////////////////////////////////////////////
// Classes
/////////////////////////////////////////////////////////////////////

struct Vector3{
	double axis[3];
	
	Vector3(){ axis[0] = 0.0; axis[1] = 0.0; axis[2] = 0.0; }
	Vector3(double x, double y, double z){ axis[0] = x; axis[1] = y; axis[2] = z; }
	void operator = (Vector3);
	void operator += (Vector3);
	void operator -= (Vector3);
	void operator *= (double);
	Vector3 operator + (Vector3);
	Vector3 operator - (Vector3);
	Vector3 operator * (double);
	double Dot(Vector3);
	Vector3 Cross(Vector3);
	double Length();
	double Distance(Vector3);
	double Normalize();
	void Dump();
};

class Efficiency{
    private:
	double *small_energy, *small_efficiency;
	double *med_energy, *med_efficiency;
	double *large_energy, *large_efficiency;
	unsigned short NsmallEff, NmedEff, NlargeEff;
	bool init_small, init_med, init_large;
	
	bool _read_eff_file(const char* fname, std::vector<double> &energy, std::vector<double> &efficiency);
	
    public:
	Efficiency();
	~Efficiency();
	
	unsigned short GetNsmall(){ return NsmallEff; }
	unsigned short GetNmedium(){ return NmedEff; }
	unsigned short GetNlarge(){ return NlargeEff; }
	
	bool IsSmallInit(){ return init_small; }
	bool IsMediumInit(){ return init_med; }
	bool IsLargeInit(){ return init_large; }
	
	unsigned short ReadSmall(const char*);
	unsigned short ReadMedium(const char*);
	unsigned short ReadLarge(const char*);
	double GetSmallEfficiency(double);
	double GetMediumEfficiency(double);
	double GetLargeEfficiency(double);
};

class AngularDist{
    private:
	std::vector<double> com_theta, dsigma_domega, integral;
	double reaction_xsection, rate;
	unsigned int num_points;
	bool init;
	
    public:
	AngularDist(){
		reaction_xsection = 0.0;
		num_points = 0;
		init = false;
	}
	~AngularDist(){
	}
	
	bool Initialize(const char*, double, double, double);
	
	double GetRate(){ return rate; }
	unsigned int GetNumPoints(){ return num_points; }
	double GetReactionXsection(){ return reaction_xsection; }
	
	double Sample();
};

class Kindeux{
    private:   	
    	double Mbeam, Mtarg, Mrecoil;
    	double Meject, Qvalue, tgt_thickness;
    	double *RecoilExStates;
    	
    	unsigned short NDist, NrecoilStates;
    	AngularDist *distributions;
    	bool ang_dist, init;
    	
    	bool get_excitations(double&);
   	void prod_proc();
    
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
    		if(ang_dist){ delete[] distributions; } 
    	}
    	
    	double GetMbeam(bool in_kg=false){ 
    		if(!in_kg){ return Mbeam; }
    		else{ return Mbeam/6.02214129E26; }
    	}
    	double GetMtarg(bool in_kg=false){ 
    		if(!in_kg){ return Mtarg; }
    		else{ return Mtarg/6.02214129E26; }
    	}
    	double GetMrecoil(bool in_kg=false){ 
    		if(!in_kg){ return Mrecoil; }
    		else{ return Mrecoil/6.02214129E26; }
    	}
    	double GetMeject(bool in_kg=false){ 
    		if(!in_kg){ return Meject; }
    		else{ return Meject/6.02214129E26; }
    	}
   	
    	void Initialize(double, double, double, double, double, unsigned short, double*, double);
    	bool SetDist(std::vector<std::string>&, double, double);
    	bool FillVars(double, double, double, double&, Vector3&);
    	bool FillVars(double, double, double, double&, double&, Vector3&, Vector3&);
	double ConvertAngle2CoM(double, double, double);
	double GetEnergies(double, double, unsigned short, double*, double*);
	void Sample(double*);
	void Print();
};

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

double min(double, double);
double max(double, double);
double frand();
short atos(const char*);
void UnitRandom(Vector3&);
double WrapValue(double, double, double);
unsigned short GetLines(const char*);
void Cart2Sphere(double, double, double, double&, double&, double&);
void Cart2Sphere(double, double, double, Vector3&);
void Cart2Sphere(Vector3, Vector3&);
double beta2(double, double);
double btoep(double);
double zeff(double, double);
double deff(double);
double shell(double);
double dedxp(double, double, double);
double dedx(double, double, double);
double range(double, double, double, double, double);
double algip(double);
void ncdedx(double, double, double, double, double, double, double&, double&, double&);
double de(double, double, double, double, double);
double linear(double, double, double, double, double);
double momentum(double, double);
double radlength(unsigned short, unsigned short);
double rndgauss0(double);
void rndgauss1(double&, double&, double&, double&, double&);
void Sphere2Cart(double, double, double, double&, double&, double&);
void Sphere2Cart(double, double, double, Vector3&);
void Sphere2Cart(Vector3, Vector3&);
double velocity(double, double);
void straggleA(double&, double, double, double, double, double);
void transform(double, double, double, double, double&, double&);
void strag_targ(double, double, double, double, double, double, double&, double&, double);
void targ_thick(double, double, double, double, double, double&);
void unitV(double, double, double, double&, double&, double&, double&);
void AngDist_read(std::string, unsigned short&, double*, double*, double&);
void SRIMread(std::string, bool&, unsigned short&, double*, double*, double*, double*, double*, bool);
std::string to_str(double);
double Interpolate(double, double, double, double, double);

#endif
