// strag_dE_cyl.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:36:26 2014

// Warning on line 7: Found subroutine strag_dE_cyl, input types cannot be determined

void strag_dE_cyl(A,Z,detThick,DetPhi,DetRad,theta_old,phi_old,energy,theta_out,phi_out,X,conv_det){ 
	
	 
	//
	// strag_dE_cyl 1.0 written by S.D.Pain on 23/02/2005
	//
	// strag_dE_cyl 1.1 modified by S.D.Pain on 7/03/2005
	// to untilise transform2.0
	//
	// Subroutine for the calcualtion of the angular straggling in
	// the dE layer of a telescope, and application of its effect
	// on the measured angles in the E detector.
	//
	//     A = A of ion
	//     Z = Z of ion
	//     DetThick = Apparent detector thickness in um (because of incident angle)
	//     DetPhi = azimuthal angle of detector
	//     theta_old = polar angle of incident ion (wrt target)
	//     phi_old = azimuthal angle of incident ion (wrt target)
	//     energy = energy of ion
	//     theta_new = polar angle of outgoing ion (wrt target)
	//     phi_new = azimuthal angle of outgoing ion (wrt target)
	//     X = radiation length of detector material
	//     conv_det = detector thickness conversion from um to mg/cm2
	//
	//     PhiDet = azimuthal angle of incidnent ion wrt detector
	//
	
	double A, Z, detThick; 
	double conv_det, DetThick_mg; 
	double pi, length1, length2, PhiDet1, PhiDet2, DetPhi; 
	double theta_old, phi_old, theta_scat, phi_scat; 
	double theta_scatW; 
	double theta_new, phi_new, energy; 
	double oldX, oldY, oldZ; 
	double newX, newY, newZ, X; 
	double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3; 
	double DetSep, detAngle, DetRad; 
	double dummy, dummylength; 
	double theta_out, phi_out; 
	double rndgauss0, rad2deg; 
	
	pi = 3.141592654d0; 
	rad2deg=180.0/pi; 
	
	// Calculate the detector thickness in mg/cm^2, and determine the
	// scattering width
	//      DetThick = 140
	//      energy = 5.0
	DetThick_mg = DetThick*conv_det; 
	
	call straggleA(theta_scatW,energy,Z,A,DetThick_mg,X); 
	
	//      print*, theta_scatW/pi*180.0
	
	theta_scat = rndgauss0(theta_scatW); 
	theta_scat = sqrt(theta_scat**2*2.0); 
	phi_scat = rand(0)*2.0*pi; 
	
	
	dummylength = 1.0; 
	// Direction of the particle
	
	// Convert the incident vector to cartesian, as this will be needed later
	call sphere2cart(dummylength,theta_old,phi_old,oldX,oldY,oldZ); 
	
	
	
	// Transform the scattered ion vector to the lab coordinates (direction relative to target)
	call transform(theta_old,phi_old,theta_scat,phi_scat,theta_new,phi_new); 
	
	
	
	// Convert the scattered vector to cartesians. This will be used to
	// calculate the angle of the scattered ion relative to the detector. (direction again)
	call sphere2cart(dummylength,theta_new,phi_new,newX,newY,newZ); 
	
	
	
	// Just to shut the compiler up!
	PhiDet2 = 0.0; 
	// Calculate the phi angle of the ejectile for the dE layer
	// (ie angle relative to plane of detector, looking along beam axis)
	// Note this calcualtion is also used in det_thick_cyl - perhaps it needs
	// to be a function?
	if ((SQRT((Phi_new-DetPhi)**2))>(0.5*pi)){
		if (DetPhi>pi){
			PhiDet2 = Phi_new-(DetPhi-2.0*pi); 
			elseif(DetPhi<=pi)then; 
			PhiDet2 = Phi_new-(DetPhi+2.0*pi); 
		} 
		} else{ 
		PhiDet2 = Phi_new-DetPhi; 
	} 
	
	DetSep = 1.0; // Separation of detectors (cm) possibly read this in at some point
	
	
	// Calculate absolute angle of the scattered ion wrt to detector plane
	// (I think this is an approximation - check)
	//      detAngle = acos(cos(PhiDet2)*cos(theta_new-(0.5*pi)))
	// I think this is the real thing
	
	
	// Calculate the length of the vector from the scattering point to the
	// point of incidence on the E detector
	
	
	// Just to shut the compiler up!
	
	
	if (DetPhi>pi){
		PhiDet1 = Phi_old-(DetPhi-2.0*pi); 
		elseif(DetPhi<=pi)then; 
		PhiDet1 = Phi_old-(DetPhi+2.0*pi); 
	} 
	} else{ 
	PhiDet1 = Phi_old-DetPhi; 
} 
 

// Calculate the length of the vector from the target to
// the scattering point
// Length due to phi

// Length due to theta



// the vector from target to scattering point
X1 = X1*length1; 
Y1 = Y1*length1; 
Z1 = Z1*length1; 



// the vector from scattering point to plane of E detector
X2 = X2*length2; 
Y2 = Y2*length2; 
Z2 = Z2*length2; 


// The vector from target to point of incidence on E detector
y3 = y1+y2; 
z3 = z1+z2; 


// Convert this back to spherical polars, and return theta and phi
& theta_out,phi_out); 



// Diagnostics
//      print*,'angle_old ',theta_old/pi*180.0,phi_old/pi*180.0
//      print*,'angle_scat ',theta_scat/pi*180.0,phi_scat/pi*180.0
//      print*,'angle_new ',theta_out/pi*180.0,phi_out/pi*180.0
//     &  ,dummylength


} 
