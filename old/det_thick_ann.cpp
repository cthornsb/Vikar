// det_thick_ann.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:37:01 2014

// Warning on line 7: Found subroutine det_thick_ann, input types cannot be determined

void det_thick_ann(Theta,Phi,DetThickness,Thickness){ 
	 
	
	//
	// det_thick_ann 1.0 written by S.D.Pain on 20/4/2006
	//
	// Subroutine to calculate the detector thickness seen by a particle
	// for annular detectors. Currently annular detectors are planar
	
	double theta, phi, DetThickness, thickness; 
	//      double precision pi
	
	//      pi = 3.141592654d0
	
	//
	// Thickness due to theta
	thickness = DetThickness/(sqrt(pow(cos(theta), 2)) ); 
	
	
	return 0; 
} 

