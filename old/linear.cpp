// linear.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:19:18 2014

double linear(double xmin, double xmax, double ymin, double ymax, double x){
	// linear 1.0 written by S.D.Pain on 24/11/2004
	// Function for linear interpolation between two points

	double grad = (ymax-ymin)/(xmax-xmin); 
	double cint = ymin - grad*xmin; 
	return x*grad + cint; 
} 

