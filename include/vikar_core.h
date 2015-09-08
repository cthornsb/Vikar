// vikar_core.h
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
extern const double rad2deg, LN2;

/////////////////////////////////////////////////////////////////////
// Classes
/////////////////////////////////////////////////////////////////////

struct Vector3{
	double axis[3];
	
	Vector3(){ axis[0] = 0.0; axis[1] = 0.0; axis[2] = 0.0; }
	Vector3(double x, double y, double z){ axis[0] = x; axis[1] = y; axis[2] = z; }
	const Vector3& operator = (const Vector3&);
	const Vector3& operator += (const Vector3&);
	const Vector3& operator -= (const Vector3&);
	const Vector3& operator *= (const double&);
	Vector3 operator + (const Vector3&) const ;
	Vector3 operator - (const Vector3&) const ;
	Vector3 operator * (const double&) const ;
	double Dot(const Vector3 &) const ;
	Vector3 Cross(const Vector3 &) const ;
	double Length() const ;
	double Distance(const Vector3 &) const ;
	double Normalize();
	std::string Dump() const ;
};

class Matrix3{
  private:
	double components[3][3];
	
	void _initialize();
	
  public:
	Matrix3();
	Matrix3(double theta_, double phi_);
	Matrix3(const Vector3 &vector_);
	void GetUnitX(Vector3 &vector){ vector.axis[0] = components[0][0];  vector.axis[1] = components[1][0];  vector.axis[2] = components[2][0]; }
	void GetUnitY(Vector3 &vector){ vector.axis[0] = components[0][1];  vector.axis[1] = components[1][1];  vector.axis[2] = components[2][1]; }
	void GetUnitZ(Vector3 &vector){ vector.axis[0] = components[0][2];  vector.axis[1] = components[1][2];  vector.axis[2] = components[2][2]; }
	Vector3 GetUnitX(){ return Vector3(components[0][0], components[1][0], components[2][0]); }
	Vector3 GetUnitY(){ return Vector3(components[0][1], components[1][1], components[2][1]); }
	Vector3 GetUnitZ(){ return Vector3(components[0][2], components[1][2], components[2][2]); }
	void SetRow1(double p1, double p2, double p3){ components[0][0] = p1; components[0][1] = p2; components[0][2] = p3; }
	void SetRow2(double p1, double p2, double p3){ components[1][0] = p1; components[1][1] = p2; components[1][2] = p3; }
	void SetRow3(double p1, double p2, double p3){ components[2][0] = p1; components[2][1] = p2; components[2][2] = p3; }
	void SetRotationMatrixSphere(double theta_, double phi_);
	void SetRotationMatrixSphere(const Vector3 &vector_);
	void SetRotationMatrixCart(double x_, double y_, double z_);
	void SetRotationMatrixCart(const Vector3 &vector_);
	void Transform(Vector3 &vector_);
	void Transpose(Vector3 &vector_);
	void Dump();
};

class AngularDist{
  private:
	double *com_theta, *dsigma_domega, *integral;
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
		if(init){
			delete[] com_theta;
			delete[] dsigma_domega;
			delete[] integral;
		}
	}
	
	bool Initialize(const char* fname, double mtarg, double tgt_thickness, double beam_current);
	bool Initialize(unsigned int num_points_, double *angle_, double *xsection_);
	
	double GetRate(){ return rate; }
	unsigned int GetNumPoints(){ return num_points; }
	double GetReactionXsection(){ return reaction_xsection; }
	
	bool Sample(double &com_angle);
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
   	
   	void Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, unsigned int NrecoilStates_, double *RecoilExStates_);
	bool SetDist(std::vector<std::string> &fnames, double total_targ_mass, double tgt_thickness_, double incident_beam_current);
	bool SetRutherford(double coeff_);
	bool FillVars(double Beam_E, double &Ejectile_E, Vector3 &Ejectile, double &comAngle, int recoil_state=-1, int solution=-1, double theta=-1);
	bool FillVars(double Beam_E, double &Ejectile_E, double &Recoil_E, Vector3 &Ejectile, Vector3 &Recoil, 
				  double &comAngle, int recoil_state=-1, int solution=-1, double theta=-1);
	double ConvertAngle2Lab(double, double, double);
	void Print();
};

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

bool IsInVector(const std::string &input_, const std::vector<std::string> &str_vector_);
void RandomCircle(double radius_, double offset_, Vector3 &beam);
void RandomGauss(double fwhm_, double offset_, Vector3 &beam);
bool SetBool(std::string input_, std::string text_, bool &output);
bool Prompt(std::string prompt_);
std::string Parse(std::string input);
void Parse(std::string input, float *arr, unsigned int num_values);
double Order(double input_);
double Dist3d(const Vector3&, const Vector3&);
double dabs(double);
double min(double, double);
double max(double, double);
double frand();
double frand(double, double);
void UnitSphereRandom(Vector3&);
void UnitSphereRandom(double&, double&);
double UnitCircleRandom();
double WrapValue(double, double, double);
unsigned int GetLines(const char*);
void Cart2Sphere(double, double, double, double&, double&, double&);
void Cart2Sphere(double, double, double, Vector3&);
void Cart2Sphere(const Vector3&, double&, double&, double&);
void Cart2Sphere(const Vector3&, Vector3&);
void Cart2Sphere(Vector3&);
double radlength(unsigned int, unsigned int);
double rndgauss0(double);
void Sphere2Cart(double, double, double, double&, double&, double&);
void Sphere2Cart(double, double, double, Vector3&);
void Sphere2Cart(const Vector3&, double&, double&, double&);
void Sphere2Cart(const Vector3&, Vector3&);
void Sphere2Cart(Vector3&);
void straggleA(double&, double, double, double, double, double);
double Interpolate(double, double, double, double, double);

#endif
