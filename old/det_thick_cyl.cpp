// det_thick_cyl.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:37:11 2014

// Warning on line 7: Found subroutine det_thick_cyl, input types cannot be determined

void det_thick_cyl(Theta,Phi,DetPhi,DetThickness,Thickness){ 
	 
	
	//
	// det_thick_cyl 1.0 written by S.D.Pain on 12/02/2005
	//
	// Subroutine to calculate the detector thickness seen by a particle
	// for detectors in cylindrical polar coordinates
	
	double theta, phi, DetPhi, DetThickness, thickness; 
	double PhiDet, pi; 
	
	pi = 3.141592654d0; 
	
	// Just to shut the compiler up!
	PhiDet = 0.0; 
	// Calculate the thickness the ejectile sees in the dE layer
	// Note this calcualtion is also used in strag_dE - perhaps it needs
	// to be a function?
	if ((SQRT((Phi-DetPhi)**2))>(0.5*pi)){
		if (DetPhi>pi){
			PhiDet = Phi-(DetPhi-2.0*pi); 
			elseif(DetPhi<pi)then; 
			PhiDet = Phi-(DetPhi+2.0*pi); 
		} 
		} else{ 
		PhiDet = Phi-DetPhi; 
	} 
	// Thickness due to phi
	thickness = DetThickness/cos(PhiDet); 
	// Thickness due to theta
	thickness = thickness/cos(theta-(0.5*pi)); 
	//	thickness = DetThickness
	
	
	return 0; 
} 

