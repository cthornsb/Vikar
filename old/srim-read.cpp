// srim-read.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:39:56 2014

#ifndef SRIM_READ_CPP
#define SRIM_READ_CPP

#include <iostream>
#include <fstream>
#include <cstring>

void SRIMread(std::string fName, bool &SRIM_stat, int Npoints, double *energy, double *dedx, double *range, double *longitude, double *latitude, bool &convert){ 
	// SRIMread 1.0 written by S.D.Pain on 24/11/2004
	//
	// Subroutine for reading in data from a SRIM output file fName
	// The material density is read in, and all data in the main table
	// are read in.
	// Conversions are made to ensure all measurements are in um
	// Ranges (and widths) are converted from um to mg/cm^2,
	// using the density read from the SRIM file
	// SRIM_stat provides limited error reporting (!)
	// SRIM_stat onwards are retunred variables	
 
	int inerror, i;  
	char *dummyC, *unit_E, *unit_range, *unit_longitude, *unit_latitude; 
	double dedxE, dedxN, density, conv;
	std::string junk;
	
	// SRIM_stat stores error status - T = good, F = bad
	SRIM_stat = true; 
	
	std::ifstream file10(fName.c_str()); 
	for(i = 1; i <= 10; i++){
		file10 >> junk; 
	} 
	
	// Read in density in g/cm^3
	file10 >> dummyC >> dummyC >> dummyC >> density; 
	density = density*1000.0; // *16.0/14.0 // Convert density to mg/cm^3
	conv = 1.0e-4*density; 
	for (i = 1; i <= 4; i++){
		file10 >> junk; 
	} 
	
	while(strcmp(dummyC, "------") != 0){
		file10 >> dummyC; 
	} 
	
	Npoints = 0; // Zero the data point counter
	
	// Read in the main data points from the SRIM output file
	while(!file10.eof() && SRIM_stat){
		Npoints=Npoints+1; 
		file10 >> energy[Npoints] >> unit_E >> dedxE >> dedxN >> range[Npoints] >> unit_range;
		file10 >> longitude[Npoints] >> unit_longitude >> latitude[Npoints] >> unit_latitude; 
		
		// Add the dedx for electric and nuclear effects
		dedx[Npoints] = dedxE + dedxN; 
		
		// Make sure the energies are in MeV
		if(strcmp(unit_E,"eV") == 0){ energy[Npoints] = energy[Npoints]/1000000.0; }
		else if(strcmp(unit_E,"keV") == 0){ energy[Npoints] = energy[Npoints]/1000.0; }
		else if(strcmp(unit_E,"MeV") == 0){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the range values are in um
		if(strcmp(unit_range, "A") == 0){ range[Npoints] = range[Npoints]*0.0001; }
		else if(strcmp(unit_range, "mm") == 0){ range[Npoints] = range[Npoints]*1000.0; }
		else if(strcmp(unit_range, "um") == 0){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the longitude values are in um
		if(strcmp(unit_longitude, "A") == 0){ longitude[Npoints] = longitude[Npoints]*0.0001; }
		else if(strcmp(unit_longitude, "mm") == 0){ longitude[Npoints] = longitude[Npoints]*1000.0; }
		else if(strcmp(unit_longitude, "um") == 0){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the latitude values are in um
		if(strcmp(unit_latitude, "A") == 0){ latitude[Npoints] = latitude[Npoints]*0.0001; }
		else if(strcmp(unit_latitude, "mm") == 0){ latitude[Npoints] = latitude[Npoints]*1000.0; }
		else if(strcmp(unit_latitude, "um") == 0){  } 
		else{ SRIM_stat = false; } 
		
		// Convert from length to mg/cm^2 if necessary
		if (convert){
			range[Npoints] = range[Npoints]*conv; 
			longitude[Npoints] = longitude[Npoints]*conv; 
			latitude[Npoints] = latitude[Npoints]*conv; 
		} 
	} // Read in data loop
	
	Npoints = Npoints-1; 
	file10.close(); 

	if(!SRIM_stat){
		std::cout << "Arse Biscuits!!";
	} 
} 

#endif
