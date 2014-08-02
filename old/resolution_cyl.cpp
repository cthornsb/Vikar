// resolution_cyl.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:42:40 2014

// Warning on line 9: Found subroutine resolution_cyl, input types cannot be determined

#include <cmath>

void resolution_cyl(energy,theta,phi,beamspot,targ_angle,Eres,Pres,PresE,length,radius,newEnergy,newTheta){ 
	
	
	 
	
	//
	// resolution 1.0 written by S.D.Pain on 20/01/2004
	//
	//
	// resolution_cyl 1.1 updated by S.D.Pain on 11/12/2013
	// to improve handling of beamspot effects
	//
	// Subroutine to smear the particles by the detector resolution
	// in energy and position
	// Pres = position resolution at 5.8MeV
	// Pres1 = position resolution at 'energy'
	
	double energy, theta, phi, Eres, Pres, Pres1, PresE; 
	double pi, radius, newEnergy, NewTheta, position; 
	double length, position1, beamspot, targ_angle; 
	double dummy_trans, dummy_long, d_radius; 
	double ran_theta; 
	
	double rndgauss0; 
	
	pi = 3.141592654d0; 
	
	// resistive ORRUBA for Kelly's 14C
	//      if(energy.gt.0.5.and.theta.le.pi/2.0)then
	//      newEnergy = energy + rndgauss0(Eres)
	//      elseif(energy.gt.1.5)then
	//      newEnergy = energy + rndgauss0(Eres)
	//      else
	//       newEnergy=0
	//      endif
	
	if (energy>0.2){
		newEnergy = energy + rndgauss0(Eres); 
		} else{ 
		newEnergy=0; 
	} 
	
	// ALL UNFINISHED - NEEDS CHECKING!!!!!
	
	//      PhiH = PhiHit - phi
	//      rad1 = radius/cos(PhiH)
	//      rad = sqrt(rad1**2 + (z - 0.5*length)**2)
	//      thetaMax = 0.5*pi - asin((z - 0.5*length)/rad)
	//      rad = sqrt(rad1**2 + (z + 0.5*length)**2)
	//      thetaMin = 0.5*pi - asin((z + 0.5*length)/rad)
	
	d_radius = rndgauss0(beamspot); // SDP beamspot size adjustment, 2013-12-11
	d_radius = abs(d_radius); // SDP beamspot size adjustment, 2013-12-11
	ran_theta = rand(0)*pi; // SDP beamspot size adjustment, 2013-12-11
	dummy_long = d_radius*std::cos(ran_theta); // SDP beamspot size adjustment, 2013-12-11
	dummy_trans = d_radius*std::sin(ran_theta); // SDP beamspot size adjustment, 2013-12-11
	
	//      dummy_trans = rndgauss0(beamspot/sqrt(2.0)) !SDP beamspot size effect adjusted 2013-12-11
	//      dummy_long  = rndgauss0(beamspot/sqrt(2.0)) !SDP beamspot size effect adjusted 2013-12-11
	
	position = (radius+dummy_trans)*tan(0.5*pi-theta); // SDP beamspot transverse size adjustment, 2013-10-15
	//      print*, theta/pi*180.0,position
	position = position + dummy_long*tan(targ_angle); // SDP beamspot longitudinal size adjustment, 2013-10-15
	
	
	
	
	Pres1 = (Pres/length*PresE)/energy*length; 
	//      Pres1 = Pres !32 !(Pres/length*PresE)/energy*length
	
	Position1 = position + rndgauss0(Pres1); // currently set to no energy-dependence!+rndgauss0(testes) !second term for beamspot size
	//      print*,''
	//      print*, position,position1
	newtheta = -(atan(Position1/radius) - 0.5*pi); 
	//      print*, theta/pi*180.0,NewTheta/pi*180.0
	
	return 0; 
} 

