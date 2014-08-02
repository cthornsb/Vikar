// resolution_plan.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:42:58 2014

#ifndef RESOLUTION_PLAN_CPP
#define RESOLUTION_PLAN_CPP

#include "vikar_lib.cpp"

void resolution_plan(double energy, double theta, double phi, double Eres, double Pres, double PresE, double DetZ_plan, bool xPlane,
		     bool yPlane, double DetXMin, double DetYMin, double width, double length, int Nstrips, int StripHit, 
		     double newEnergy, double newTheta, double newPhi){ 
	// resolution_plan (VANDLE) 1.0 adapted by S.D.Pain on 2014/01/22
	// Subroutine to smear the particles by the detector resolution
	// in energy and position, for planar detectors
	
	// Pres = position resolution at 5.8MeV
	// Pres1 = position resolution at 'energy'
	
	double Pres1; 
	double pi; 
	double dummy; 
	double hit_r, hit_x, hit_y; 
	
	double t_res; 
	double vel, n_mass, time, hit_vect; 
	
	pi = 3.1415926540; 
	
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

#endif
