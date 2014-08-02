// de.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 20:30:46 2014

double de(double dx, double emass, double epart, double zpart, double sp){ 
	// to calculate energy lost over dx cms assuming a quadratic relationship between e & x.
	double deltae = sp*dx; 
	double enew = epart-deltae; 
	
	if(enew <= 0.0){ return epart; }
	
	spges = dedx(emass, enew, zpart); 
	return dx*(0.750*sp+(0.250*spges/sp)*spges); 
} 
