// AngDist_read.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:32:45 2014

#include <fstream>

void AngDist_read(const char* fName, int &Npoints, double *angle, double *integral, double &max_integral){ 
	// AngDist_read 1.0 written by S.D.Pain on 5/05/2006
	//
	// Subroutine for reading in an angular distribution profile from
	// from fName.
	// An angular distribution profile is a cumulative integration of
	// the angular distribution.
	// The input file should be of the form:
	// [angle (deg)] [cumulative integral(0-angle)]
	// where the angles must span the range 0 to 180 degrees.	
	
	int i; 
	max_integral = 0.0; 
	// DetSet_stat stores error status - T = good, F = bad

	std::ifstream file10(fName); 
	Npoints = 0; // Zero the points counter
	
	// Read in the main data points from the SRIM output file
	while(!file10.eof()){
		Npoints = Npoints+1; 
		file10 >> angle[Npoints] >> integral[Npoints]; 
	} // Read in data loop
	
	if(Npoints > 0){ Npoints = Npoints-1; }
	file10.close(); 
	
	for (i = 0; i < Npoints; i++){
		max_integral = integral[Npoints]; 
	} 
} 

#endif
