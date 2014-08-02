// strag_targ.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:41:11 2014

void strag_targ(A,Z,targ_thick,theta_old,phi_old,energy,theta_new,phi_new,X){ 
	// strag_targ 1.0 written by S.D.Pain on 20/01/2004
	//
	// strag_targ 1.1 modified by S.D.Pain on 7/03/2005
	// to untilise transform2.0
	//
	// Subroutine to calculate the angular straggling of an ion in
	// the target. The A,Z of the ion are read in, along with the
	// theta,phi and energy of the ion. The average radiation length
	// of the target material is read in as X.
	//
	// The calculation of the width of the scattering distribution
	// is calculated by straggleA. A Gaussian weighted scattering angle
	// based on this width is calculated, using rndgauss0
	//
	// The new theta,phi to which the ion is scattered is returned.
	
	double A, Z, targ_thick; 
	double pi, theta_scatW, X; 
	double theta_old, phi_old, theta_scat, phi_scat; 
	double theta_new, phi_new, energy; 
	double rndgauss0; 
	
	pi = 3.1415926540; 
	
	// Calculate the straggling width
	straggleA(theta_scatW,energy,Z,A,targ_thick,X); 
	
	// Select the scattering angle of the ion wrt its initial direction
	theta_scat = rndgauss0(theta_scatW); 
	theta_scat = std::sqrt(pow(theta_scat, 2)*2.0); 
	phi_scat = rand(0)*2.0*pi; 
	
	// Determine the absolute lab angle to which the ion is scattered
	transform(theta_old,phi_old,theta_scat,phi_scat,theta_new,phi_new); 
} 
