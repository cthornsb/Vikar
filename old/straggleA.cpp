// straggleA.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:40:49 2014

void straggleA(theta,energy,Z,A,thickness,X){ 
	// straggleA 1.0 written by S.D.Pain on 20/01/2004
	//
	// Subroutine to calculate the width of a gaussian distribution
	// of angles from the straggling of an energetic ion in a medium
	//   theta = sigma of distribution (spatial)
	//   Energy = energy of particle
	//   thickness = thickness of material
	//   X = radiation length of stopping material
	//   A = Mass number of ion
	//   Z = charge of ion
	
	// CURRENTLY ONLY TESTED FOR A LIMITED RANGE OF IONS, ENERGIES and TARGETS
	double theta, v, p, A, Z, thickness, X; 
	double energy, velocity, momentum; 
	
	// Dignostics - to be removed whan satisfied
	//      energy = 5.80
	//      Z = 1.0 !4
	//      A = 1.0
	//      thickness = 0.2 !25.000
	//
	//      v = sqrt(2.0*energy/(A*931.5))
	//      p = sqrt(2.0*A*931.5*energy)
	
	// Most of this (the v and p calcualtion) cancels out in the equation
	// It's left in for clarity, but I might remove it, as it's pointless
	v = velocity(energy,A); 
	p = momentum(energy,A); 
	
	// Dignostic - to be removed whan satisfied
	theta = 13.6/(v*p)*z*std::sqrt(thickness/X)*(1.0+0.038*std::log(thickness/X)); 
	theta = theta*std::sqrt(2.0); 
} 


