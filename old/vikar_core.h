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

#include "planar.h"
//#include "cylindrical.h"
//#include "annular.h"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

extern const double pi, deg2rad, rad2deg;

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

double min(double, double);

double max(double, double);

double frand();

short atos(const char*);

double WrapValue(double, double, double);

unsigned short GetLines(const char*);

void cart2sphere(double, double, double, double&, double&, double&);

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

void sphere2cart(double, double, double, double&, double&, double&);

double velocity(double, double);

void straggleA(double&, double, double, double, double, double);

void transform(double, double, double, double, double&, double&);

void strag_targ(double, double, double, double, double, double, double&, double&, double);

void targ_thick(double, double, double, double, double, double&);

void unitV(double, double, double, double&, double&, double&, double&);

void AngDist_read(std::string, unsigned short&, double*, double*, double&);

void SRIMread(std::string, bool&, unsigned short&, double*, double*, double*, double*, double*, bool);

//void kindeux(double, double, double, double, double, double, double, double, double, double, unsigned short, unsigned short, double, double*, double*, double&, double&, double&, double&, double&, double&);

//void product_proc(double, double, double, double, double, double, double, double, double, double*, double*, double, unsigned short, unsigned short*, double**, double**, double**, double**, double*, double*, double*, double*, double*, unsigned short, unsigned short*, double*, double**, double**, double**, double**, unsigned short, unsigned short*, double**, double**, double**, double**, double*, double*, double, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, bool*, bool*, double*, double*, double*, double*, double*, double*, double*, double*, unsigned short&, unsigned short&, unsigned short&, unsigned short&, unsigned short&, unsigned short&, unsigned short&, unsigned short&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, bool&);

class Kindeux{
    private:   	
    	double Mbeam, Mtarg, Mrecoil, Meject;
    	
    	bool ang_dist;
    	unsigned short NDistPoints;
    	double Int_max;
    	double *DistAng, *DistInt;
    
    public:
    	Kindeux(){
    		ang_dist = false;
    		Mbeam = 0.0; Mtarg = 0.0;
    		Mrecoil = 0.0; Meject = 0.0;
    	}
    	~Kindeux(){
    		if(ang_dist){
    			delete[] DistAng;
    			delete[] DistInt;
    		}
    	}
    	
    	void Initialize(double, double, double, double);
    	void SetDist(unsigned short, double, double*, double*);    	
    	void FillVars(double, double, double, double, double, double, double&, double&, double&, double&, double&, double&);
};

class ProdProc{
    private:
	double A, Z, beamspot, thickness;
	double targ_depth, targ_angle, targRadL;
	double *E_targ, *Erange_targ;

	unsigned short Ndet_plan, *Nstrips_plan;
	double *DetZ_plan, *E_det, *Erange_det, *dEthick_plan, *Ethick_plan;
	double **DetXMax, **DetXMin, **DetYMax, **DetYMin; 
	double conv_det, DetRadL;
	
	double *Eres_dE_plan, *Pres_dE_plan, *Eres_E_plan, *Pres_E_plan; 
	double *Pres_dE_en_plan, *Pres_E_en_plan, *detWidth_plan, *detLength_plan;
	bool *xPlane, *yPlane; 
	bool init;

    public:
    	ProdProc(){
    		init = false;
    	}
    	~ProdProc(){
    	}
    	
    	void Initialize(double, double, double, double, double, double, double*, double*, double, unsigned short, 
    			unsigned short*, double*, double**, double**, double**, double**, double*, double*, double, double, double*, double*, 
    			double*, double*, double*, double*, double*, double*, bool*, bool*, double*, double*);
    	bool FillVars(double,double, double, unsigned short&, unsigned short&, double&, double&, double&, double&, 
    		      double&, double&, double&, double&, double&, double&, double&, double&, double&, bool&);
};

#endif
