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
extern const double rad2deg, LN2, avagadro;

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

struct Matrix3{
	double components[3][3];
	
	Matrix3();
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
	
	unsigned int NDist, NrecoilStates;
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
   	
	void Initialize(double, double, double, double, double, unsigned int, double*, double);
	bool SetDist(std::vector<std::string>&, double, double);
	bool FillVars(double, double&, Vector3&);
	bool FillVars(double, double&, double&, Vector3&, Vector3&);
	bool FillVars(double, double, double*, double*, double*, double*, double*, double*, double*);
	double ConvertAngle2Lab(double, double, double);
	double GetEnergies(double, double, unsigned int, double*, double*);
	void Sample(double*);
	void Print();
};

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

std::string to_str(double);
std::string Parse(std::string input);
void Parse(std::string input, float *arr, unsigned int num_values);
double Dist3d(const Vector3&, const Vector3&);
double dabs(double);
double min(double, double);
double max(double, double);
double frand();
void UnitSphereRandom(Vector3&);
void UnitSphereRandom(double&, double&);
double UnitCircleRandom();
double WrapValue(double, double, double);
unsigned int GetLines(const char*);
void Cart2Sphere(double, double, double, double&, double&, double&);
void Cart2Sphere(double, double, double, Vector3&);
void Cart2Sphere(const Vector3&, Vector3&);
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
double momentum(double, double);
double radlength(unsigned int, unsigned int);
double rndgauss0(double);
void rndgauss1(double&, double&, double&, double&, double&);
void Sphere2Cart(double, double, double, double&, double&, double&);
void Sphere2Cart(double, double, double, Vector3&);
void Sphere2Cart(const Vector3&, Vector3&);
double velocity(double, double);
void straggleA(double&, double, double, double, double, double);
void transform(double, double, double, double, double, double);
double Interpolate(double, double, double, double, double);

#endif
