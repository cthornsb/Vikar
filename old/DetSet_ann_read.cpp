// DetSet_ann_read.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:58 2014

// Warning on line 9: Found subroutine DetSet_ann_read, input types cannot be determined

#include <fstream>

void DetSet_ann_read(fName,Ndet,Nstrips,DetInner,DetOuter,DetdPhi,DetPhi,DetZ,dEthick,Ethick,Eres_dE,Eres_E){ 
	
	
	 
	//
	// DetSet_ann_read 1.0 written by S.D.Pain on 7/4/2005
	//
	// Subroutine for reading in annular detector details and locations
	// from fName.
	// The annular detector coordinates are specified in XXXXX?
	// and are read in in cm from the origin.
	// The input file format should be of the form:
	//
	// (Detector Inner Radius (cm)) (Detector Outer Radius (cm))
	// (Detector Phi Coverage (deg)) (Detector Phi Centre (deg))
	// (N. Strips)  (Z position (cm))
	// (dE thickness (um))   (E thickness (um))
	// (dE Energy Resolution (MeV))  (E Energy Resolution (MeV))
	
	int Ndet, Nstrips[200]; 
	double DetInner[200], DetOuter[200]; 
	double DetPhi[200], DetdPhi[200], DetZ[200]; 
	double dEthick[200], Ethick[200]; 
	double Eres_dE[200], Eres_E[200]; 
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
		file10 >> detInner(Ndet) >> DetOuter(Ndet) >> DetdPhi(Ndet) >> DetPhi[Ndet] >> Nstrips[Ndet] >> DetZ(Ndet) >> dEthick(Ndet) >> Ethick(Ndet) >> Eres_dE(Ndet) >> Eres_E(Ndet); 
		
		
		
		
		DetdPhi(Ndet) = DetdPhi(Ndet)*deg2rad; 
		DetPhi[Ndet] = DetPhi[Ndet]*deg2rad; 
	} // Read in data loop
	
	
	if (Ndet>0){Ndet=Ndet-1; }
	file10.close(); 
	
	
	return 0; 
} 
