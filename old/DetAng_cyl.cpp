// DetAng_cyl.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:00 2014

// Warning on line 8: Found subroutine DetAng_cyl, input types cannot be determined
// Warning on line 36: Cannot determine order of arguments of do loop

void DetAng_cyl(length,width,Nstrips,radius,z,phi,ThetaMax,ThetaMin,PhiMax,PhiMin){ 
	
	
	 
	//
	// DetAng_cyl 1.0 written by S.D.Pain on 9/01/2005
	//
	// Subroutine for calculating the range of theta and phi (spherical
	// polar coordinates) spanned by each strip of a strip detector
	// oriented in cylindrical polar coordinates, with the strip direction
	// parallel to the z axis.
	// The length, width and strip division are passed from the parent program,
	// along with the radius, z location and angle (phi) of the detector, in
	// cylindrical polars.
	// The maximum and minumum theta and phi angles for each strip are calculated
	// and returned.
	
	double length, width, radius, z, phi; 
	double rad, rad1, rad2, rad3; 
	double PhiMax[50], PhiMin[50], xmin, xmax; 
	double ThetaMax[50], thetaMin[50]; 
	
	double pi; 
	int Nstrips, n; 
	
	pi = 3.141592654d0; 
	
	
	for (n = 1; n <= Nstrips; n++){
		// Calculate the position (cm) of the edges of the strip
		xmax = (width/Nstrips)*n - width/2.0; 
		xmin = (width/Nstrips)*(n-1) - width/2.0; 
		
		// rad1,2 = radii to the strip edge, in the xy plane
		rad1 = SQRT(pow((xmax), 2 )+pow( radius, 2) ); 
		PhiMax[n] = phi + asin(xmax/rad1); 
		
		rad2 = SQRT(pow((xmin), 2 )+pow( radius, 2) ); 
		PhiMin(n) = phi + asin(xmin/rad2); 
		
		// The usual bloody conversions to keep 0<phi<2*pi
		if (PhiMax[n]<0){ PhiMax[n] = 2.0*pi + PhiMax[n]; }
		if (PhiMin(n)<0){ PhiMin(n) = 2.0*pi + PhiMin(n); }
		
		if (PhiMax[n]>(2.0*pi)){ PhiMax[n] = PhiMax[n] - 2.0*pi; }
		if (PhiMin(n)>(2.0*pi)){ PhiMin(n) = PhiMin(n) - 2.0*pi; }
		
		//        print*, 'phi min= ',phiMin(n)/3.14159*180.0
		//        print*, 'phi max= ',phiMax(n)/3.14159*180.0
		
		// Find the strip edge nearest the detector centre
		// This will be the closest to the target, so has the greatest
		// coverage in theta
		if (rad1<rad2){
			rad3 = rad1; 
			} else{ 
			rad3 = rad2; 
		} 
		
		// Calculate the radius to the strip corner to get thetaMax
		rad = sqrt(pow(rad3, 2 )+pow( (z - 0.5*length), 2) ); 
		thetaMax(n) = 0.5*pi - asin((z - 0.5*length)/rad); 
		//        print*, 'theta max= ',thetaMax(n)/3.14159*180.0
		
		// And for ThetaMin
		rad = sqrt(pow(rad3, 2 )+pow( (z + 0.5*length), 2) ); 
		thetaMin[n] = 0.5*pi - asin((z + 0.5*length)/rad); 
		//        print*, 'theta min= ',thetaMin(n)/3.14159*180.0
		//        PhiMax(n) = phi + asin(xmax/rad1)
		
		
		//	thetaMin(n) = (z - 0.5*length)
	} 
	
	return 0; 
} 

