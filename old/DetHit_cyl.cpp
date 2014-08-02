// DetHit_cyl.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:43:14 2014

// Warning on line 7: Found subroutine DetHit_cyl, input types cannot be determined

void DetHit_cyl(length,width,radius,z,phi,thetaHit,PhiHit,Hit){ 
	
	
	 
	//
	// DetHit_cyl 1.0 written by S.D.Pain on 9/01/2005
	//
	// Subroutine for fine calculation of hits.
	// The precise range of theta (spherical polar coordinates)
	// spanned by the detector at the phi of the hit is calculated
	// to verify the hit.
	// The length, width and strip division are passed from the parent program,
	// along with the radius, z location and angle (phi) of the detector, in
	// cylindrical polars.
	
	
	double length, width, radius, z, phi; 
	double rad, rad1, thetaMax, ThetaMin; 
	double ThetaHit, PhiHit, PhiH; 
	double pi; 
	int hit; 
	
	pi = 3.141592654d0; 
	
	PhiH = PhiHit - phi; 
	rad1 = radius/cos(PhiH); 
	rad = sqrt(pow(rad1, 2 )+pow( (z - 0.5*length), 2) ); 
	thetaMax = 0.5*pi - asin((z - 0.5*length)/rad); 
	rad = sqrt(pow(rad1, 2 )+pow( (z + 0.5*length), 2) ); 
	thetaMin = 0.5*pi - asin((z + 0.5*length)/rad); 
	
	if (ThetaHit<=thetaMax && ThetaHit>=thetaMin){
		
		hit = hit; 
		} else{ 
		hit = hit - 1; 
	} 
	
	return 0; 
} 

