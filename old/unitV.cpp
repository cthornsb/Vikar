// unitV.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:42:15 2014

void unitV(xl,yl,zl,x,y,z,length){ 
	// unitV 1.0 written by S.D.Pain on 20/11/2004
	//
	// Subroutine to read in a cartesian vector (xl,yl,zl)
	// and calculate and return its unit vector (x,y,z) and length
	
	double length, xl, yl, zl, x, y, z, norm; 
	
	// Put in check of length > 0    ?????
	length = std::sqrt(pow(xl, 2)+pow(yl, 2)+pow(zl, 2)); 
	norm = std::sqrt(1.0/(pow(xl, 2)+pow(yl, 2)+pow(zl, 2))); 
	x = xl*norm; 
	y = yl*norm; 
	z = zl*norm; 
} 
