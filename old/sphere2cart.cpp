// sphere2cart.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:39:42 2014

#include <cmath>

void sphere2cart(r,theta,phi,x,y,z){ 
	// sphere2cart 1.0 written by S.D.Pain on 4/12/2004
	//
	// Subroutine for converting a vector from spherical polar coordinates
	// to cartesian coordinates
	//
	// (r,theata,phi) are passed in, and (x,y,z) are calculated and
	// returned
	
	double theta, phi, r, pi; 
	double x, y, z, temp; 
	
	pi = 3.1415926540; 
	
	z = r*std::cos(theta); 
	temp = r*std::sin(theta); 
	x = temp*std::sin(phi); 
	y = temp*std::cos(phi); 
	
	return 0; 
} 

