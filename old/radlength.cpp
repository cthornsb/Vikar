// radlength.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:40:46 2014

double functionradlength(int A, int Z){
	// radlength 1.0 written by S.D.Pain on 11/02/2005
	// Function to calculate the radiation length of a material
	// in mg/cm^2
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135

	return 7.164d5*A/(Z*(Z+1.0)*std::log(287.0/(std::sqrt(Z)))); 
} 


