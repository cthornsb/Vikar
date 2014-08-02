// DetSet_cyl_read.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:36:25 2014

// Warning on line 9: Found subroutine DetSet_cyl_read, input types cannot be determined

#include <fstream>

void DetSet_cyl_read(fName,Ndet,Nstrips,detLength,DetWidth,DetRad,DetZ,DetPhi,dEthick,Ethick,Eres_dE,Pres_dE,Eres_E,Pres_E,Pres_dE_energy,Pres_E_energy){ 
	
	
	
	 
	//
	// DetSet_read 1.0 written by S.D.Pain on 9/01/2005
	//
	// Subroutine for reading in detector details and locations
	// from fName.
	// The detector coordinates are specified in cylindrical polars,
	// and are read in in degrees.
	// The input file format should be of the form:
	//
	// (Detector length (cm)) (Detector Width (cm))  (N. Strips)  (Radius (cm))
	// (Z position (cm))   (Phi (deg))   (dE thickness (um))   (E thickness (um))
	// (dE Energy Resolution (MeV))  (dE Pos Resolution (cm))
	// (E Energy Resolution (MeV))  (E Pos Resolution (cm))
	// (dE Pos Resolution Energy (MeV))  (E Pos Resolution Energy (MeV))
	
	int Ndet, Nstrips[200]; 
	double detLength[200], DetWidth[200]; 
	double DetRad[200], DetZ[200], DetPhi[200]; 
	double dEthick[200], Ethick[200]; 
	double Eres_dE[200], Pres_dE[200], Eres_E[200]; 
	double Pres_E[200], Pres_dE_energy[200]; 
	double Pres_E_energy[200]; 
	double deg2rad; 
	
	double pi; 
	int inerror; 
	char* fName; 
	
	pi = 3.141592654d0; 
	deg2rad = pi/180.0; 
	
	// DetSet_stat stores error status - T = good, F = bad
	
	std::ifstream file10(fName); 
	
	
	Ndet = 0; // Zero the detector counter
	
	// Read in the main data points from the SRIM output file
	while (inerror==0){
		
		Ndet=Ndet+1; 
		file10 >> detLength(Ndet) >> DetWidth(Ndet) >> Nstrips[Ndet] >> DetRad(Ndet) >> DetZ[Ndet] >> DetPhi(Ndet) >> dEthick(Ndet) >> Ethick(Ndet) >> Eres_dE(Ndet) >> Pres_dE(Ndet) >> Eres_E(Ndet) >> Pres_E(Ndet) >> Pres_dE_energy(Ndet) >> Pres_E_energy(Ndet); 
		
		
		
		
		
		DetPhi(Ndet) = DetPhi(Ndet)*deg2rad; 
		
	} // Read in data loop
	
	
	if (NDet>0){ Ndet=Ndet-1; }
	file10.close(); 
	
	
	return 0; 
} 
