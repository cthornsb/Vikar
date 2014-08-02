// DetSet_plan_read.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:36:39 2014

#include <fstream>

void DetSet_plan_read(fName,Ndet,Nstrips){ 
	// DetSet_plan_read 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for reading in planar detector details and locations
	// from fName.
	// The planar detector coordinates are specified in cartesians,
	// and are read in in cm from the origin.
	// The input file format should be of the form:
	//
	// (Detector length (cm)) (Detector Width (cm))  (N. Strips)  (Z position (cm))
	// (X Position (cm)) (Y Position (cm))  (plane of strips ('x' or 'y')   (dE thickness (um))   (E thickness (um))
	// (dE Energy Resolution (MeV))  (dE Pos Resolution (cm))
	// (E Energy Resolution (MeV))  (E Pos Resolution (cm))
	// (dE Pos Resolution Energy (MeV))  (E Pos Resolution Energy (MeV))
	
	// returns detLength, DetWidth, DetZ, DetX, DetY, xPlane, yPlane, dEthick, Ethick, Eres_dE, Pres_dE, Eres_E, Pres_E, Pres_dE_energy, Pres_E_energy
	
	int Ndet, Nstrips[200]; 
	double detLength[200], DetWidth[200]; 
	double DetZ[200], DetX[200], DetY[200]; 
	double dEthick[200], Ethick[200]; 
	double Eres_dE[200], Pres_dE[200], Eres_E[200]; 
	double Pres_E[200]; 
	double Pres_dE_energy[200], Pres_E_energy[200]; 
	double deg2rad; 
	
	double pi; 
	int inerror; 
	bool xPlane[200], yPlane[200]; 
	char* fName; 
	character*1 plane; 
	
	pi = 3.141592654d0; 
	deg2rad = pi/180.0; 
	
	// DetSet_stat stores error status - T = good, F = bad
	
	std::ifstream file10(fName); 
	
	Ndet = 0; // Zero the detector counter
	
	// Read in the main data points from the SRIM output file
	while(!file10.eof()){
		Ndet=Ndet+1; 
		file10 >> detLength[Ndet] >> DetWidth[Ndet] >> Nstrips[Ndet] >> DetZ[Ndet];
		file10 >> DetX[Ndet] >> DetY[Ndet] >> plane >> dEthick[Ndet] >> Ethick[Ndet];
		file10 >> Eres_dE[Ndet] >> Pres_dE[Ndet] >> Eres_E[Ndet] >> Pres_E[Ndet];
		file10 >> Pres_dE_energy[Ndet] >> Pres_E_energy[Ndet]; 

		if (plane == "x"){
			xPlane[Ndet] = true; 
			yPlane[Ndet] = false; 
		else if(plane == "y"){
			xPlane[Ndet] = false; 
			yPlane[Ndet] = true; 
		} 	
	} // Read in data loop
	
	if(Ndet > 0){ Ndet=Ndet-1; }
	
	file10.close(); 
} 
