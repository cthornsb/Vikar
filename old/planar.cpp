// Det_Planar.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:20 2014

#include "planar.h"
#include "vikar_core.h"

/////////////////////////////////////////////////////////////////////
// DetAng_plan.f
/////////////////////////////////////////////////////////////////////

void DetAng_plan(double length, double width, unsigned short Nstrips, double z, double x, double y, 
		 bool xplane, bool yplane, double *xMax, double *xMin, double *yMax, double *yMin){ 
	// DetAng_plan 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for calculating the range of x and y spanned by each strip
	// of a strip detector planar coordinates (normal parallel to the z axis),
	// with the strip direction parallel to the x or y axis.
	// The length, width and strip division are passed from the parent program,
	// along with the z location of the plane, the xy coordinates of the detector
	// centre, and the direction of the strips.
	// The maximum and minumum x and y coordinates for each strip are calculated
	// and returned.
	
	// "---------------" << std::endl;
	// "Detector: " << detNo << " Strip: " << n << std::endl;
	// "DetXMax =" << "\t" << DetXMax[detNo][n] << std::endl;
	// "DetXMin =" << "\t" << DetXMin[detNo][n] << std::endl;
	// "DetYMax =" << "\t" << DetYMax[detNo][n] << std::endl;
	// "DetYMin =" << "\t" << DetYMin[detNo][n] << std::endl;
	
	unsigned short n;  
	for(n = 0; n < Nstrips; n++){
		// Calculate the position (cm) of the edges of the strip
		if(xplane){
			yMax[n] = (width/Nstrips)*n - width/2.0 + y; 
			yMin[n] = (width/Nstrips)*(n-1) - width/2.0 + y; 
			xMax[n] = x+(length/2.0); 
			xMin[n] = x-(length/2.0); 
		}
		else if(yplane){
			xMax[n] = (width/Nstrips)*n - width/2.0 + x; 
			xMin[n] = (width/Nstrips)*(n-1) - width/2.0 + x; 
			yMax[n] = y+(length/2.0); 
			yMin[n] = y-(length/2.0); 
		} 
	} 
} 

/////////////////////////////////////////////////////////////////////
// DetHit_plan.f
/////////////////////////////////////////////////////////////////////

void DetHit_plan(double **DetXMax, double **DetXMin, double **DetYMax, double **DetYMin, unsigned short Ndet_plan, unsigned short *Nstrips_plan,
		 double *DetZ_plan, double theta, double phi, unsigned short &hit, unsigned short &DetHit, unsigned short &StripHit){ 
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
	 
	double hit_r, hit_x, hit_y; 
	unsigned short n, s;  	 
	
	for(n = 0; n < Ndet_plan; n++){
		for(s = 0; s < Nstrips_plan[n]; s++){
			// For each detector, check the particle is heading towards it
			if (((theta*rad2deg) < 90.0 && DetZ_plan[n] > 0) || ((theta*rad2deg) > 90.0 && DetZ_plan[n] < 0)){
				hit_r = DetZ_plan[n]*std::tan(theta); // radius of hit from z-axis in detector plane
				hit_x = hit_r*std::sin(phi); // x position of hit in detector plane
				hit_y = hit_r*std::cos(phi); // y position of hit in detector plane
		
				// Check if particle is within detector area
				if (hit_x <= DetXMax[n][s] && hit_x >= DetXMin[n][s] && hit_y <= DetYMax[n][s] && hit_y >= DetYMin[n][s]){
					hit = hit + 1; 
					StripHit = s; 
					DetHit = n; 
				} 
			} // theta check
		} // Strip
	} // Detector
}

/////////////////////////////////////////////////////////////////////
// DetSet_plan_read.f
/////////////////////////////////////////////////////////////////////

void DetSet_plan_read(const char* fName, unsigned short Ndet, unsigned short *Nstrips, double *detLength, double *detWidth, double *DetZ, double *DetX, double *DetY,
		      bool *xPlane, bool *yPlane, double *dEthick, double *Ethick, double *Eres_dE, double *Pres_dE, double *Eres_E,
		      double *Pres_E, double *Pres_dE_energy, double *Pres_E_energy){ 
	// DetSet_plan_read 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for reading in planar detector details and locations
	// from fName.
	// The planar detector coordinates are specified in cartesians,
	// and are read in in cm from the origin.
	// The input file format should be of the form:
	//
	// (Detector length (cm)) (Detector Width (cm))  (N. Strips)  (Z position (cm))
	// (X Position (cm)) (Y Position (cm))  (plane of strips ('x' or 'y')   (dE thickness (um))   (E thickness (um))
	// (dE Energy Resolution (MeV))  (dE Pos Resolution (cm))
	// (E Energy Resolution (MeV))  (E Pos Resolution (cm))
	// (dE Pos Resolution Energy (MeV))  (E Pos Resolution Energy (MeV))
	std::string plane;

	// DetSet_stat stores error status - T = good, F = bad
	std::ifstream file10(fName); 	
	
	// Read in the main data points from the SRIM output file
	for(unsigned short i = 0; i < Ndet; i++){
		file10 >> detLength[i] >> detWidth[i] >> Nstrips[i] >> DetZ[i];
		file10 >> DetX[i] >> DetY[i] >> plane >> dEthick[i] >> Ethick[i];
		file10 >> Eres_dE[i] >> Pres_dE[i] >> Eres_E[i] >> Pres_E[i];
		file10 >> Pres_dE_energy[i] >> Pres_E_energy[i]; 

		if(plane == "x"){
			xPlane[i] = true; 
			yPlane[i] = false;
		} 
		else if(plane == "y"){
			xPlane[i] = false; 
			yPlane[i] = true; 
		} 	
	} // Read in data loop
	
	file10.close(); 
} 

/////////////////////////////////////////////////////////////////////
// det_thick_plan.f
/////////////////////////////////////////////////////////////////////

void det_thick_plan(double theta, double phi, double DetThickness, double &thickness){ 
	// det_thick_plan 1.0 written by S.D.Pain on 8/6/2005
	//
	// Subroutine to calculate the detector thickness seen by a particle
	// for detectors in planar coordinates
	
	// Thickness due to theta
	thickness = DetThickness/(std::sqrt(pow(std::cos(theta), 2.0)) ); 
}

/////////////////////////////////////////////////////////////////////
// resolution_plan.f
/////////////////////////////////////////////////////////////////////

void resolution_plan(double energy, double theta, double phi, double Eres, double Pres, double PresE, double DetZ_plan, bool xPlane,
		     bool yPlane, double DetXMin, double DetYMin, double width, double length, unsigned short Nstrips, unsigned short StripHit, 
		     double newEnergy, double newTheta, double newPhi){ 
	// resolution_plan (VANDLE) 1.0 adapted by S.D.Pain on 2014/01/22
	// Subroutine to smear the particles by the detector resolution
	// in energy and position, for planar detectors
	
	// Pres = position resolution at 5.8MeV
	// Pres1 = position resolution at 'energy'
	
	double Pres1; 
	double dummy; 
	double hit_r, hit_x, hit_y; 
	
	double t_res; 
	double vel, n_mass, time, hit_vect; 
	
	// **** begin VANDLE specific ******
	n_mass = 1.0; 
	t_res = 3.0; 						// beam timing resolution in ns (Changed to 3 ns for vandle timing by Cory)
	vel = velocity(energy,n_mass); 
	// **** end VANDLE specific ******
	
	hit_r = DetZ_plan*tan(theta); 				// radius of hit from z-axis in detector plane
	hit_x = hit_r*sin(phi); 				// x position of hit in detector plane
	hit_y = hit_r*cos(phi); 				// y position of hit in detector plane
	
	// **** begin VANDLE specific ******
	hit_vect = std::sqrt(pow(DetZ_plan, 2)+pow(hit_r, 2)); 
	time = (((frand()-0.5)*5.0)+hit_vect)/vel; 		// 5.0 = bar thickness in cm (large bar)
	// time = (((rand(0)-0.5)*3.0)+hit_vect)/vel 		// 3.0 = bar thickness in cm (medium & small bar)
	time = time +rndgauss0(t_res); 
	// time = hit_r/vel 					// 2.5 = bar thickness in cm
	vel = hit_vect/time; 
	// **** end VANDLE specific ******
	
	// newEnergy = energy + rndgauss0(Eres) 		// removed for VANDLE
	newEnergy = 0.5*n_mass*pow(vel, 2 ); 			// added for VANDLE
	
	// Pres1 = (Pres/length*PresE)/energy*length 		// Pres at 'energy' !removed for VANDLE
	Pres1 = Pres;						// added for VANDLE
	
	if(xPlane){
		// hit_y = detyMin + (width/Nstrips*StripHit)-(width/Nstrips)/2.0	//assumes detyMin is bottom of detector, not strip
		hit_x = hit_x + rndgauss0(Pres1);
		hit_y = DetYMin + (width/Nstrips)/2.0; 
	}
	else if(yPlane){ 
		// hit_x = detxMin + (width/Nstrips*StripHit)-(width/Nstrips)/2.0 	// as above
		hit_x = DetXMin + (width/Nstrips)/2.0; 
		hit_y = hit_y + rndgauss0(Pres1); 
	} 
	
	cart2sphere(hit_x,hit_y,DetZ_plan,dummy,newTheta,newPhi); 
} 

/////////////////////////////////////////////////////////////////////
// strag_dE_plan.f
/////////////////////////////////////////////////////////////////////

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

	// Separation of detectors (cm) possibly read this in at some pounsigned short 
	DetSep = 1.0;

	// Calculate absolute angle of the scattered ion wrt to detector plane 
	if(theta_new <= pi/2.0){ detAngle = theta_new; }
	else if(theta_new > pi/2.0){ detAngle = pi-theta_new; }

	// Calculate the length of the vector from the scattering pounsigned short to the
	// pounsigned short of incidence on the E detector
	length2 = DetSep/cos(detAngle);

	// Calculate the length of the vector from the target to
	// the scattering point
	length1 = DetZ/cos(theta_old);

	// the vector from target to scattering point
	unitV(oldX,oldY,oldZ,X1,Y1,Z1,dummy);
	X1 = X1*length1;
	Y1 = Y1*length1;
	Z1 = Z1*length1;

	// the vector from scattering pounsigned short to plane of E detector
	unitV(newX,newY,newZ,X2,Y2,Z2,dummy);
	X2 = X2*length2;
	Y2 = Y2*length2;
	Z2 = Z2*length2;

	// The vector from target to pounsigned short of incidence on E detector
	X3 = X1+X2;
	Y3 = Y1+Y2;
	Z3 = Z1+Z2;
	 
	// Convert this back to spherical polars, and return theta and phi
	cart2sphere(X3,Y3,Z3,dummylength,theta_out,phi_out);
}
