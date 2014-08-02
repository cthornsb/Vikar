// resolution_ann.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:39:16 2014

// Warning on line 8: Found subroutine resolution_ann, input types cannot be determined


void resolution_ann(energy,theta,phi,Eres,DetZ,DetInner,DetOuter,DetPhi,Nstrips,StripHit,newEnergy,newTheta,newPhi){ 
	
	
	
	 
	
	
	//
	// resolution_ann 1.0 written by S.D.Pain on 20/04/2006
	//
	// Subroutine to smear the particles by the detector resolution
	// in energy and position, for annular detectors
	
	
	// Pres = position resolution at 5.8MeV
	// Pres1 = position resolution at 'energy'
	
	double energy, theta, phi, Eres, detPhi; 
	double DetZ, DetInner, DetOuter; 
	double pi, newEnergy, NewTheta, NewPhi; 
	double hit_r, r_width; 
	
	double rndgauss0; 
	int Nstrips, StripHit; 
	
	pi = 3.141592654d0; 
	
	if (energy>0.2){
		newEnergy = energy + rndgauss0(Eres); 
		} else{ 
		newEnergy=0.0; 
	} 
	
	
	// Radius (in cylindrical polars) to the centre of the strip
	hit_r = DetInner+((DetOuter-DetInner)/Nstrips)*(float(striphit)-0.5); 
	
	r_width = ((DetOuter-DetInner)/Nstrips); 
	
	//      print*, hit_r, r_width
	
	if (DetZ>0){
		NewTheta = datan((hit_r+(rand(0)-0.5)*r_width)/DetZ); 
		elseif(DetZ<0)then; 
		NewTheta = pi+datan((hit_r+(rand(0)-0.5)*r_width)/DetZ); 
	} 
	
	newPhi = detPhi; 
	//      print*, theta,newTheta
	
	//       newtheta = theta
	
	return 0; 
} 

