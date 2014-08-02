// transform.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:41:42 2014

void transform(theta1,phi1,theta2,phi2,theta,phi){ 
	// transform 2.0 written by S.D.Pain on 4/03/2005
	//
	// Subroutine for transforming the a spherical polar vector
	// from one refernce frame to another.
	// (theta2,phi2) is a vector in the master frame
	// (theta1,phi1) is measured relative to (theta2,phi2).
	// (theta,phi) is (theta1,phi1) in the master frame
	
	double theta1, phi1, theta2, phi2, theta, phi; 
	double term1, term2, temp, x1, y1, x2, y2, x, y; 
	double beamX, beamY, beamZ, dummy, pi; 
	double dumtheta1, dumphi1, dumtheta2, dumphi2; 
	bool swap; 
	
	swap = false; 
	
	pi = 3.1415926540; 
	dummy = 1.0; 
	
	// copy the input angles to different variables, and use the copies in
	// the subroutine, as they get modified.
	dumtheta1 = theta1; 
	dumphi1 = phi1; 
	dumtheta2 = theta2; 
	dumphi2 = phi2; 
	
	// Check whether the vector is pointing backward of 90 degrees (polar)
	// If so, reflect its direction around, so that it points forwards.
	// The transformation can then be computed, and the vector reflected
	// back again. This avoids edge-of-the-world effects.
	
	if (dumtheta1 > (0.5*pi)){
		swap = true; 
		
		// Doesn't appear that any transformation is needed here - perhaps
		// worth checking, though...
		//        call sphere2cart(dummy,theta2,phi2,beamX,beamY,beamZ)
		//        beamX = -beamX
		//        beamY = -beamY
		//        beamZ = -beamZ
		//        call cart2sphere(beamX,beamY,beamZ,dummy,theta2,phi2)
		
		sphere2cart(dummy,dumtheta1,dumphi1,beamX,beamY,beamZ); 
		beamX = -beamX; 
		beamY = -beamY; 
		beamZ = -beamZ; 
		cart2sphere(beamX,beamY,beamZ,dummy,dumtheta1,dumphi1); 
	} 
	
	// Calculate the total angle between the two vectors. This is the
	// effective polar angle.
	term1 = dumtheta1 + dumtheta2*(cos(dumphi2-dumphi1)); 
	term2 = dumtheta2*(sin((dumphi2-dumphi1))); 
	theta = sqrt(pow(term1, 2 )+pow( term2, 2) ); 
	
	x1 = dumtheta1*sin(dumphi1); 
	y1 = dumtheta1*cos(dumphi1); 
	x2 = dumtheta2*sin(dumphi2); 
	y2 = dumtheta2*cos(dumphi2); 
	
	x = x1+x2; 
	y = y1+y2; 
	
	temp = x/(sqrt(pow(x, 2)+pow(y, 2)) ); 
	phi = asin(temp); 
	
	if (x>=0.0){ phi = acos(y/(std::sqrt(pow(x, 2)+pow(y, 2)))); } 
	else{ phi = 2.0*3.14159-acos(y/(sqrt(pow(x, 2)+pow(y, 2)))); }
	
	// Keeps theta & phi within limits. Not needed, with reflection
	// procedure
	//      if(theta.gt.3.14159)then
	//       theta = 2.0*3.14159-theta
	//        phi = phi+3.14159
	//        if(phi.gt.2.0*3.14159) phi = phi-2.0*3.14159
	//      endif
	
	// If a reflection was made, reflect back again.
	if (swap){
		sphere2cart(dummy,theta,phi,beamX,beamY,beamZ); 
		beamX = -beamX; 
		beamY = -beamY; 
		beamZ = -beamZ; 
		cart2sphere(beamX,beamY,beamZ,dummy,theta,phi); 
	} 
	
	return 0; 
} 
