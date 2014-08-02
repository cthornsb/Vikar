// kindeux.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:37:47 2014

#ifndef KINDEUX_CPP
#define KINDEUX_CPP

#include <cmath>

#include "vikar_lib.cpp"

void kindeux(double Ebeam, double theta_beam, double phi_beam, double Q, double ExEject, double ExRecoil, double Mbeam, double Mtarg, double Mrecoil, double Meject, int Dist_type,
	     int NDistPoints, double Int_max, double *DistAng, double *DistInt, double &theta_Recoil, double &phi_Recoil, double &theta_Eject, double &phi_Eject, double &Erecoil, double &Eeject){ 
	// kindeux 1.0 written by S.D.Pain on 24/11/2004
	//
	// kindeux 2.0 updated by S.D.Pain on 4/3/2005
	//    - Updated to utilise transform3
	//    - Cleaned up old junk, diagnostics etc
	//
	// kindeux 3.0 updated by S.D.Pain on 5/5/2006
	//    - Updated to employ non-isotropic angular distributions
	//
	// Subroutine for calculating two body reaction kinematics
	// The reaction energy is read in, along with the direction of the incident
	// ion. The Q value, excitation energies and particle masses are read in, also.
	// The reaction occurs evenly distributed over all solid angle, in the CoM frame,
	// or is defined by an angular distribution profile for the ejectile, passed
	// from the parent program. This profile is in the form of a cumulative integral
	// of an angular distribution of the form dSigma/dTheta (not dSigma/dOmega) in
	// the CoM system, with 0 degrees defined by normal kinematics.
	// The lab energies and directions are calculated for the recoil and ejectile.
	// The recoil and ejectile directions are altered due to the angle of the
	// incoming beam particle. E, theta, phi are returned for recoil and ejectile.
	
	// Variables passed in 
	double temp; 
	int i; 
	
	// Internal Variables
	double E, vBeam, vCoM, vRecoil, vEject; 
	double recoilX, recoilY, recoilZ; 
	double ejectX, ejectY, ejectZ; 
	double theta_beam_com, phi_beam_com, lab_theta_beam; 
	double lab_Phi_beam; 
	double beamX, beamY, beamZ; // lab direction of the beam
	double CoMX, CoMY, CoMZ; 
	double deg2rad; 

	deg2rad=pi/180.0; 
	
	// Calculate beam velocity
	vBeam = velocity(Ebeam,Mbeam);
	 
	// Calculate the CoM velocity in the lab
	vCoM = vBeam*Mbeam/(Mbeam+Mtarg); 
	
	// Set the direction of the CoM frame to along the Z axis.
	// Perform all calcualtions in this frame, and rotate the
	// lab vectors later.
	theta_beam_com = 0.0; 
	phi_beam_com = 0.0; 
	
	// Store the direction of the beam particle in the lab
	// to be used later to determine the real angles of the
	// recois and ejectile
	lab_theta_beam = theta_beam; 
	lab_Phi_beam = phi_beam; 
		
	// Convert the CoM velocity to cartesions (should be entirely along
	// the Z direction), and convert to a unit vector. Proabably
	// pointless now, as everything is aligned along the Z axis
	sphere2cart(vCoM,theta_beam_com,phi_beam_com,beamX,beamY,beamZ); 
	unitV(beamX,beamY,beamZ,CoMX,CoMY,CoMZ,vCoM); 
	
	// Randomly select outgoing coordinates in CoM frame
	// I think this samples evenly over solid angle... check
	// for the ejectile
	// NB theta = 0 corresponds to forward angles in the lab
	// so for the inverse kinematics case, need to reverse the
	// coordinates so theta = 0 gives backward angles in the lab
	phi_Eject = 2.0*pi*frand(); 
	
	if(Dist_type == 1){
		i = 1; 
		temp = frand()*Int_max; 
		while(temp > DistInt[i]){
			i=i+1; 
		} 
		theta_Eject = linear(DistInt[i-1],DistInt[i],DistAng[i-1],DistAng[i],temp); 
		
		theta_Eject = theta_Eject*deg2rad; 
		// swap directions for inverse kinematics cases
		if (Mbeam>Mtarg){ theta_Eject = pi-theta_Eject; }
		} else{ 
		theta_Eject = acos(-2.0*frand()+1.0); 
	} 
	
	// Calcualte the recoil's theta and phi
	theta_Recoil = pi-theta_Eject; 
	if (phi_Eject<pi){ phi_Recoil = phi_Eject+pi; }
	if (phi_Eject>=pi){ phi_Recoil = phi_Eject-pi; }
	
	// Calculate the CoM energy, and add on the Q value and subtract
	// any excitations
	E = 0.5*Mbeam*pow((vBeam-vCoM), 2 )+ 0.5*Mtarg*pow(vCoM, 2 ); 
	E = E + Q - ExEject - ExRecoil; 
	
	// Calculate the magnitudes of the velocities of the reaction products
	// in the CoM frame
	vRecoil = std::sqrt((2.0*E) /(Mrecoil+pow(Mrecoil, 2)/Meject)); 
	vEject = std::sqrt((2.0*E) /(Meject+pow(Meject, 2)/Mrecoil)); 
	
	// Convert the recoil's and ejectile's CoM velocity vectors from
	// spherical polars to cartesians, and transform into the laboratory frame
	sphere2cart(vRecoil,theta_Recoil,phi_Recoil,recoilX,recoilY,recoilZ); 
	sphere2cart(vEject,theta_Eject,phi_Eject,ejectX,ejectY,ejectZ); 
	
	recoilX = recoilX + vCoM*CoMX; 
	recoilY = recoilY + vCoM*CoMY; 
	recoilZ = recoilZ + vCoM*CoMZ; 
	ejectX = ejectX + vCoM*CoMX; 
	ejectY = ejectY + vCoM*CoMY; 
	ejectZ = ejectZ + vCoM*CoMZ; 
	
	// Convert the lab frame cartesians to lab frame spherical polars
	// NB now theta_* and phi_* contain lab angles, not CoM angles
	// These variables are returned to the Master
	cart2sphere(ejectX,ejectY,ejectZ,vEject,theta_Eject,phi_Eject); 
	cart2sphere(recoilX,recoilY,recoilZ,vRecoil,theta_Recoil,phi_Recoil); 
	
	// Calculate the laboratory energy
	Erecoil = 0.5*Mrecoil*pow(vRecoil, 2.0); 
	Eeject = 0.5*Meject*pow(vEject, 2.0); 
	
	// Rotate the velocity vectors due to the incident angle of the beam particle
	transform(lab_theta_beam,lab_Phi_beam,theta_Eject,phi_Eject,theta_Eject,phi_Eject); 	
	transform(lab_theta_beam,lab_Phi_beam,theta_Recoil,phi_Recoil,theta_Recoil,phi_Recoil); 
} 

#endif
