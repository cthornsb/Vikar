// DetAng_plan.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:20 2014

void DetAng_plan(double length, double width, int Nstrips, double z, double x, double y, bool xplane,bool yplane){ 
	// DetAng_plan 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for calculating the range of x and y spanned by each strip
	// of a strip detector planar coordinates (normal parallel to the z axis),
	// with the strip direction parallel to the x or y axis.
	// The length, width and strip division are passed from the parent program,
	// along with the z location of the plane, the xy coordinates of the detector
	// centre, and the direction of the strips.
	// The maximum and minumum x and y coordinates for each strip are calculated
	// and returned.
	
	// Returns xmax, xmin, ymax, ymin
	
	double length, width, z, x, y; 
	double xMax[50], xMin[50]; 
	double yMax[50], yMin[50]; 
	
	double pi; 
	int Nstrips, n; 
	bool xplane, yplane; 
	
	pi = 3.141592654d0; 
	
	for(n = 0; n < Nstrips; n++){
		// Calculate the position (cm) of the edges of the strip
		if(xplane){
			ymax[n] = (width/Nstrips)*n - width/2.0 + y; 
			ymin[n] = (width/Nstrips)*(n-1) - width/2.0 + y; 
			xmax[n] = x+(length/2.0); 
			xmin[n] = x-(length/2.0); 
		}
		else if(yplane){
			xmax[n] = (width/Nstrips)*n - width/2.0 + x; 
			xmin[n] = (width/Nstrips)*(n-1) - width/2.0 + x; 
			ymax[n] = y+(length/2.0); 
			ymin[n] = y-(length/2.0); 
		} 
	} 
	
	return 0; 
} 

