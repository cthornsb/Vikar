// strag_dE_plan.cpp
// Feb. 14th 2014 14:38

#ifndef STRAG_DE_PLAN_CPP
#define STRAG_DE_PLAN_CPP

#include "vikar_lib.cpp"

void strag_dE_plan(double A, double Z, double detThick, double DetZ, double theta_old, double phi_old, double energy, 
		   double &theta_out, double &phi_out, double X, double conv_det){
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
	//     detThick = Apparent detector thickness in um (because of incident angle)
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

	double detThick_mg;
	double length1,length2;
	double theta_scat,phi_scat;
	double theta_scatW;
	double theta_new,phi_new;
	double oldX,oldY,oldZ;
	double newX,newY,newZ;
	double X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3;
	double DetSep,detAngle;
	double dummy,dummylength;

	// Calculate the detector thickness in mg/cm^2, and determine the
	// scattering width
	// detThick = 140 
	// energy = 5.0
	detThick_mg = detThick*conv_det;
	straggleA(theta_scatW,energy,Z,A,detThick_mg,X);
	theta_scat = rndgauss0(theta_scatW);
	theta_scat = sqrt(theta_scat*theta_scat*2.0);
	phi_scat = frand()*2.0*pi;
	dummylength = 1.0;

	// Direction of the particle
	// Convert the incident vector to cartesian, as this will be needed later
	sphere2cart(dummylength,theta_old,phi_old,oldX,oldY,oldZ);    

	// Transform the scattered ion vector to the lab coordinates (direction relative to target)    
	transform(theta_old,phi_old,theta_scat,phi_scat,theta_new,phi_new);
	     
	// Convert the scattered vector to cartesians. This will be used to
	// calculate the angle of the scattered ion relative to the detector. (direction again)
	sphere2cart(dummylength,theta_new,phi_new,newX,newY,newZ);

	// Separation of detectors (cm) possibly read this in at some point 
	DetSep = 1.0;

	// Calculate absolute angle of the scattered ion wrt to detector plane 
	if(theta_new <= pi/2.0){ detAngle = theta_new; }
	else if(theta_new > pi/2.0){ detAngle = pi-theta_new; }

	// Calculate the length of the vector from the scattering point to the
	// point of incidence on the E detector
	length2 = DetSep/cos(detAngle);

	// Calculate the length of the vector from the target to
	// the scattering point
	length1 = DetZ/cos(theta_old);

	// the vector from target to scattering point
	unitV(oldX,oldY,oldZ,X1,Y1,Z1,dummy);
	X1 = X1*length1;
	Y1 = Y1*length1;
	Z1 = Z1*length1;

	// the vector from scattering point to plane of E detector
	unitV(newX,newY,newZ,X2,Y2,Z2,dummy);
	X2 = X2*length2;
	Y2 = Y2*length2;
	Z2 = Z2*length2;

	// The vector from target to point of incidence on E detector
	X3 = X1+X2;
	Y3 = Y1+Y2;
	Z3 = Z1+Z2;
	 
	// Convert this back to spherical polars, and return theta and phi
	cart2sphere(X3,Y3,Z3,dummylength,theta_out,phi_out);
}

#endif
