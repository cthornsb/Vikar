// DetHit_plan.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:14:58 2014

void DetHit_plan(DetXMax,DetXMin,DetYMax,DetYMin){ 
	// DetHit_plan 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for determination of hits in planar detectors.
	// The precise range in x and y (cartesian coordinates)
	// spanned by the detector at the detector plane is calculated
	// to verify the hit.
	// The DetXMax,DetXMin,DetYMax,DetYMin are passed from the parent program,
	// along with the Ndet_plan, Nstrips_plan, z location and angle (theta,phi)
	// of the particle.
	// Hit, DetHit and StripHit are returned
	
	// hit_x, hit_y returned (CORY)
	
	double DetXMax[200][50], DetXMin[200][50]; 
	double DetYMax[200][50], DetYMin[200][50]; 
	double DetZ_plan[200]; 
	double Theta, Phi; 
	double hit_r, hit_x, hit_y; 
	double xMax[50], xMin[50]; 
	double yMax[50], yMin[50]; 
	double pi; 
	int hit, Nstrips, n, s, StripHit, DetHit; 
	int Nstrips_plan[200], Ndet_plan; 	
	
	pi = 3.1415926540; 
	
	for(n = 0; n < Ndet_plan; n++){
		for(s = 0; s < Nstrips_plan[n]; s++){
			// For each detector, check the particle is heading towards it
			if (((theta/pi*180.0) < 90.0 && DetZ_plan[n] > 0)||((theta/pi*180.0) > 90.0 && DetZ_plan[n] < 0)){
				hit_r = DetZ_plan[n]*tan(theta); // radius of hit from z-axis in detector plane
				hit_x = hit_r*sin(phi); // x position of hit in detector plane
				hit_y = hit_r*cos(phi); // y position of hit in detector plane
		
				// Check if particle is within detector area
				if (hit_x <= DetXMax[n][s] && hit_x >= DetxMin[n][s] && hit_y <= DetyMax[n][s] && hit_y >= DetyMin[n][s]){
					hit = hit + 1; 
					StripHit = s; 
					DetHit = n; 
				} 
			} // theta check
		} // Strip
	} // Detector
} 

