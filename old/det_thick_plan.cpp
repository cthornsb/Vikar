// det_thick_plan.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:37:21 2014

void det_thick_plan(Theta,Phi,DetThickness,Thickness){ 
	// det_thick_plan 1.0 written by S.D.Pain on 8/6/2005
	//
	// Subroutine to calculate the detector thickness seen by a particle
	// for detectors in planar coordinates
	
	double theta, phi, DetThickness, thickness; 

	// Thickness due to theta
	thickness = DetThickness/(std::sqrt(pow(std::cos(theta), 2.0)) ); 
} 

