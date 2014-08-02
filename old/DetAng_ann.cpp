// DetAng_ann.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:34:33 2014

// Warning on line 8: Found subroutine DetAng_ann, input types cannot be determined
// Warning on line 28: Cannot determine order of arguments of do loop

void DetAng_ann(DetInner,DetOuter,DetdPhi,DetPhi,Nstrips,DetZ,ThetaMax,ThetaMin,PhiMax,PhiMin){ 
	
	
	 
	//
	// DetAng_ann 1.0 written by S.D.Pain on 7/4/2006
	//
	// Subroutine for calculating the range of theta and phi (spherical
	// polar coordinates) spanned by each strip of an annular strip detector.
	
	double DetInner, DetOuter, DetdPhi, DetPhi, DetZ; 
	double PhiMax[50], PhiMin[50], rmin, rmax; 
	double ThetaMax[50], ThetaMin[50]; 
	
	double pi; 
	int Nstrips, n; 
	
	pi = 3.141592654d0; 
	
	
	for (n = 1; n <= Nstrips; n++){
		// Calculate the radius (cm) in cylindrical polars of the edges of the strip
		rmax = DetInner+((DetOuter-DetInner)/Nstrips)*n; 
		rmin = DetInner+((DetOuter-DetInner)/Nstrips)*(n-1); 
		
		// Determine the range of theta coverage for the detector, allowing for
		// forward and backward hemisphere location
		if (DetZ>0){
			ThetaMin[n] = datan(rmin/DetZ); 
			ThetaMax(n) = datan(rmax/DetZ); 
			elseif(DetZ<0)then; 
			ThetaMin[n] = pi+datan(rmax/DetZ); 
			ThetaMax(n) = pi+datan(rmin/DetZ); 
		} 
		
		// Calculate the Phi coverage of the detector
		PhiMax[n] = DetPhi + 0.5*DetdPhi; 
		PhiMin(n) = DetPhi - 0.5*DetdPhi; 
		
		// The usual bloody conversions to keep 0<phi<2*pi
		if (PhiMax[n]<0){ PhiMax[n] = 2.0*pi + PhiMax[n]; }
		if (PhiMin(n)<0){ PhiMin(n) = 2.0*pi + PhiMin(n); }
		
		if (PhiMax[n]>(2.0*pi)){ PhiMax[n] = PhiMax[n] - 2.0*pi; }
		if (PhiMin(n)>(2.0*pi)){ PhiMin(n) = PhiMin(n) - 2.0*pi; }
		
		//        print*, 'phi min= ',phiMin(n)/3.14159*180.0
		//        print*, 'phi max= ',phiMax(n)/3.14159*180.0
		
	} 
	
	return 0; 
} 

