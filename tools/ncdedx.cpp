//test.cpp

#include <iostream>
#include <cmath>

// Constants for use in stopping power calcuations
const double pi = 3.14159;
const double electron_RME = 0.510998928; // MeV
const double proton_RME = 938.272046; // MeV
const double neutron_RME = 939.565378; // MeV
const double bohr_e_radius = 2.817940326836615e-15; // m
const double e = 1.602176565E-19; // C
const double h = 4.135667516e-21; // MeV * s
const double coeff = 4*pi*pow(bohr_e_radius, 2.0)*electron_RME; // m^2 * MeV
const double alpha = h*h*9E16*bohr_e_radius/pi; // MeV^2 * m^3
const double avagadro = 6.0221413e+23; // 1/mol

// Constant coefficients for the shell correction term
const double a1 = 0.422377E-6, a2 = 3.858019E-9; // eV
const double b1 = 0.030403E-6, b2 = -0.1667989E-9; // eV
const double c1 = -0.00038106E-6, c2 = 0.00157955E-9; // eV

// Ionization potentials for H(eV) He    Li    Be    B     C     N     O     F      Ne     Na     Mg     Al
const double ionpot[13] = {18.7, 42.0, 39.0, 60.0, 68.0, 78.0, 99.5, 98.5, 117.0, 140.0, 150.0, 157.0, 163.0};

class DEDX{
  private:
	double *N, *A, *Z, *I; // Target material constants
	double average_Z, average_A;
	
	double calculate_avgz(); // Calculate the average atomic charge (unitless) for the material
	double calculate_avgA(); // Calculate the average atomic mass (amu) for the material
	
  public:
	DEDX();
	DEDX(unsigned int num_elements_);
	DEDX(unsigned int num_elements_, double *N_, double *A_, double *Z_);
	
	bool Set(unsigned int num_elements_, double *N_, double *A_, double *Z_);
	bool SetElement(unsigned int element_, double N_, double A_, double Z_);
	
	// Return the effective atomic charge (unitless) for a given value of beta_ and Z_
	double Zeff(double beta_, double Z_){ return Z_*(1-std::exp(-125.0*beta_/pow(Z_, 2.0/3.0))); }

	// Return the value of beta^2 (unitless) for a particle with a given energy_ and mass_
	// energy_ in MeV and mass_ in MeV/c^2
	double Beta2(double energy_, double mass_){ return std::sqrt(2*energy_/mass_); }

	// Return the electron density (1/m^3) for the material
	// material_density_ in g/cm^3
	double edens(double material_density_){ return avagadro*GetAverageZ()*(material_density_*(1E6))/GetAverageA(); }

	// Return the ionization potential (MeV) for a particle with a given atomic charge Z_
	// Z_ is the atomic number of the ion of interest
	double GetIonPot(unsigned int Z_);

	// Return the natural log of the average ionization potential for the material
	double lnIbar();

	// Return the shell correction term for a particle with a given energy_ in this material
	// energy_ in MeV
	double ShellCorrect(double energy_);

	// Return the density effect correction term (unitless) for a particle with a given energy_ in this material
	// energy_ in MeV
	double DensityEffect(double energy_);

	// Return the stopping power for a proton with a given energy_ in this material
	// energy_ in MeV
	double Pstop(double energy_);

	// Return the stopping power (MeV/m) for a particle with given energy_, Z_, and mass_ in this material
	// energy_ in MeV and mass_ in MeV/c^2
	double StopPower(double energy_, double Z_, double mass_);

	// Return the range (m) for a particle with given energy_, Z_, and mass_, in this material
	// energy_ in MeV and mass_ in MeV/c^2
	double Range(double energy_, double Z_, double mass_);
};

// Return the average atomic charge (unitless) for the material
double DEDX::calculate_avgZ(){
	double numerator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		numerator += N[i]*Z[i];
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		denominator += N[i];
	}
	return average_Z = numerator/denominator;
}

// Return the average atomic mass (amu) for the material
double DEDX::calculate_avgA(){
	double numerator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		numerator += N[i]*A[i];
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		denominator += N[i];
	}
	return average_A = numerator/denominator;
}

// Return the ionization potential (MeV) for a particle with a given atomic charge Z_
// Z_ is the atomic number of the ion of interest
double DEDX::GetIonPot(unsigned int Z_){ 
	if(Z_ <= 13){ return ionpot[Z_-1]*(1E-6); }
	return((9.76*Z_ + 58.8*pow(1.0*Z_, -0.19))*(1E-6));
}

// Return the natural log of the average ionization potential for the material
double DEDX::lnIbar(){
	double numerator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		numerator += N[i]*Z[i]*I[i];
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < 2; i++){
		denominator += N[i]*Z[i];
	}
	return numerator/denominator;
}

// Return the shell correction term for a particle with a given energy_ in this material
// energy_ in MeV
double DEDX::ShellCorrect(double energy_){
	if(energy_ >= 8.0){
		double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
		double f1 = ((a1/nu2) + (b1/(nu2*nu2)) + (c1*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		double f2 = ((a2/nu2) + (b2/(nu2*nu2)) + (c2*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		double Ibar = std::exp(lnIbar()); 
		return ((f1*pow(Ibar, 2.0) + f2*pow(Ibar, 3.0))/GetAverageZ());
	}
	else{
		// Aluminum and water stopping power
		double beta2 = Beta2(energy_, 1.0);
		double Zbar = GetAverageZ();
		double x = 18769.0*beta2/Zbar; // 137^2 = 187690
		double output = std::log(2*electron_RME*beta2) - lnIbar();
		
		// Need to calculate L(x)
		if(Zbar >= 13.0){ // Use L2(x) (proton stopping in Al)
		}
		else if(Zbar <= 3.33){ // Use L1(x) (proton stopping in H20)
		}
		else{ // Interpolate between L1 and L2
		}
		return 0.0;
	}
	return -1;
}

// Return the density effect correction term (unitless) for a particle with a given energy_ in this material
// energy_ in MeV
double DEDX::DensityEffect(double energy_){
	double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
	double output = std::log(alpha) + std::log(edens(1.063)) + 2.0*(std::log(nu2) - lnIbar()) - 1;
	if(output < 0.0){ return 0.0; }
	return output;
}

// Return the stopping power for a proton with a given energy_ in this material
// energy_ in MeV
double DEDX::Pstop(double energy_){
	double beta2 = Beta2(energy_, proton_RME);
	double beta = std::sqrt(beta2);
	double output = edens(1.063)*coeff*pow(Zeff(beta, 1.0), 2.0)/beta2; // C^4/(MeV*cm^3)
	output *= -1.0*(std::log(2.0*electron_RME*(beta2/(1-beta2)))-beta2-lnIbar()-ShellCorrect(energy_)-0.5*DensityEffect(energy_)); // 
	return output;
}

// Return the stopping power (MeV/m) for a particle with given energy_, Z_, and mass_ in this material
// energy_ in MeV and mass_ in MeV/c^2
double DEDX::StopPower(double energy_, double Z_, double mass_){
	double beta = std::sqrt(Beta2(energy_, mass_));
	double output = pow((Zeff(beta, Z_)/Zeff(beta, 1.0)), 2.0);
	output *= Pstop(energy_*(proton_RME/mass_));
	return output;
}

// Return the range (m) for a particle with given energy_, Z_, and mass_, in this material
// energy_ in MeV and mass_ in MeV/c^2
double DEDX::Range(double energy_, double Z_, double mass_){
	double Eci = 0.5*mass_*0.04*pow(Z_, 2.0/3.0);
	double sum = 0.0;
	if(energy_ <= Eci){
		double e1, e2;
		double step = energy_/1000.0;
		for(unsigned int i = 1; i < 1000; i++){
			e1 = step*i; e2 = step*(i+1);
			sum += 0.5*(1.0/StopPower(e2, Z_, mass_)+1.0/StopPower(e1, Z_, mass_))*(e2-e1);
		}
	}
	else{
		double e1, e2;
		double step = Eci/1000.0;
		
		// Ri(Eci, Zi, mi)
		for(unsigned int i = 1; i < 1000; i++){
			e1 = step*i; e2 = step*(i+1);
			sum += 0.5*(1.0/StopPower(e2, Z_, mass_)+1.0/StopPower(e1, Z_, mass_))*(e2-e1);
		}
		
		// Rp(Ei*mp/mi)
		double sum2 = 0.0;
		step = (energy_*proton_RME/mass_);
		for(unsigned int i = 1; i < 1000; i++){
			e1 = step*i; e2 = step*(i+1);
			sum2 += 0.5*(1.0/Pstop(e2)+1.0/Pstop(e1))*(e2-e1);
		}
		
		// Rp(Eci*mp/mi)
		double sum3 = 0.0;
		step = (Eci*proton_RME/mass_);
		for(unsigned int i = 1; i < 1000; i++){
			e1 = step*i; e2 = step*(i+1);
			sum3 += 0.5*(1.0/Pstop(e2)+1.0/Pstop(e1))*(e2-e1);
		}
		
		sum += (mass_/(proton_RME*Z_*Z_))*(sum2 - sum3);
	}
	return sum;
}

int main(){
	/*double estep = 20.0/1000.0;
	for(unsigned int i = 0; i <= 1000; i++){
		std::cout << estep*i << "\t" << Range(estep*i, 2.0, 2*proton_RME+2*neutron_RME) << std::endl;
	}*/
	std::cout << " coeff: " << coeff << std::endl;
	std::cout << " alpha: " << alpha << std::endl;
	std::cout << " Beta2: " << Beta2(30.3, (3*neutron_RME+4*proton_RME)) << std::endl;
	std::cout << " Zeff: " << Zeff(30.3, 4.0) << std::endl;
	std::cout << " GetAverageZ: " << GetAverageZ() << std::endl;
	std::cout << " edens: " << edens(1.063) << std::endl;
	std::cout << " Get Ion Pot: " << GetIonPot(4.0) << std::endl;
	std::cout << " lnIbar: " << lnIbar() << std::endl;
	std::cout << " Density Effect: " << DensityEffect(30.3) << std::endl;
	std::cout << " Shell Correct: " << ShellCorrect(30.3) << std::endl;
	std::cout << " Pstop: " << Pstop(30.3) << std::endl;
	std::cout << " StopPower: " << StopPower(30.3, 4.0, (3*neutron_RME+4*proton_RME)) << std::endl;
	//std::cout << " Range: " << Range(1.0, 1.0, proton_RME) << std::endl;

	return 0;
}
