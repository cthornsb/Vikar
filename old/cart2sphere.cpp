// cart2sphere.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:33:21 2014

#include <cmath>

void cart2sphere(x,y,z,r,theta,phi){ 
	// cart2sphere 1.0 written by S.D.Pain on 4/12/2004
	//
	// Subroutine for converting a vector from cartesian coordinates
	// to spherical polar coordinates
	//
	// (x,y,z) are passed in, and (r,theata,phi) are calculated and
	// returned
	
	double theta, phi, r, pi; 
	double x, y, z, temp; 
	
	temp = 0.0; 
	pi = 3.1415926540; 
	
	// Check the vector has length
	if(x==0.0 && y==0.0 && z==0.0){
		r = 0.0; 
		theta = 0.0; 
		phi = 0.0; 
	}
	else{ // if so, calculate (r,theta,phi)
		r = std::sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2)); 
		theta = std::acos(z/r); 
		if(x==0.0 && y==0.0){ phi = 0.0; }
		else{ temp = std::sqrt(pow(x, 2)+pow(y, 2)); } 
		
		if(x>=0.0){ phi = dacos(y/temp); }
		if(x<0.0){ phi = 2.0*pi -dacos(y/temp); } // So 0<phi<2*pi
	} 
	
	return 0; 
} 

