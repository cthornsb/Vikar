// materials.cpp
// Cory Thornsberry
// Wed Oct 8, 2014

#include "vikar_core.h"
#include "materials.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals 
/////////////////////////////////////////////////////////////////////

// Constants for use in stopping power calcuations
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

/////////////////////////////////////////////////////////////////////
// RangeTable
/////////////////////////////////////////////////////////////////////

bool RangeTable::_initialize(unsigned int num_entries_){
	if(use_table){ return false; }
	num_entries = num_entries_;
	energy = new double[num_entries];
	range = new double[num_entries];
	use_table = true;
	return true;
}

RangeTable::RangeTable(unsigned int num_entries_){
	_initialize(num_entries_);
}

RangeTable::~RangeTable(){
	if(use_table){
		delete[] energy;
		delete[] range;
	}
}
	
bool RangeTable::Init(unsigned int num_entries_){
	return _initialize(num_entries_);
}

bool RangeTable::Init(unsigned int num_entries_, double startE, double stopE, double A, double Z, Material *mat){
	if(!_initialize(num_entries_)){ return false; }
	
	// Use Material to fill the arrays
	double step = (stopE-startE)/(num_entries-1);
	double mass = (Z*proton_RME+(A-Z)*neutron_RME);
	for(unsigned int i = 0; i < num_entries; i++){
		energy[i] = startE + i*step; // Energy in MeV
		range[i] = mat->Range(energy[i], Z, mass); // Range in m
		//std::cout << energy[i] << "\t" << range[i] << std::endl;
	}

	return true;
}

bool RangeTable::Set(unsigned int pt, double energy_, double range_){
	if(!use_table){ return false; }
	if(pt > num_entries){ return false; }
	energy[pt] = energy_;
	range[pt] = range_;
	return true;
}

double RangeTable::GetRange(double energy_){
	if(!use_table){ return -1; }
	for(unsigned int i = 0; i < num_entries-1; i++){
		if(energy_ == energy[i]){ return range[i]; }
		else if(energy_ == energy[i+1]){ return range[i+1]; }
		else if(energy_ >= energy[i] && energy_ <= energy[i+1]){
			// Interpolate and return the result
			return (((range[i+1]-range[i])/(energy[i+1]-energy[i]))*(energy_-energy[i])+range[i]);
		}
	}
	return -1;
}

double RangeTable::GetEnergy(double range_){
	if(!use_table){ return -1; }
	for(unsigned int i = 0; i < num_entries-1; i++){
		if(range_ == range[i]){ return energy[i]; }
		else if(range_ == range[i+1]){ return energy[i+1]; }
		else if(range_ >= range[i] && range_ <= range[i+1]){
			// Interpolate and return the result
			return (((energy[i+1]-energy[i])/(range[i+1]-range[i]))*(range_-range[i])+energy[i]);
		}
	}
	return -1;
}

/////////////////////////////////////////////////////////////////////
// Efficiency
/////////////////////////////////////////////////////////////////////

// Read NewVIKAR VANDLE efficiency file and load data into arrays
// Returns the number of data points which are found in the file
// Assumes the following efficiency file format for each point in file
// Lab_Energy(MeV) Efficiency
bool Efficiency::_read_eff_file(const char* fname, std::vector<double> &energy, std::vector<double> &efficiency){
	std::ifstream eff_file(fname);
	if(!eff_file.good()){ return false; }
	float values[2];

	while(true){
		for(unsigned int i = 0; i < 2; i++){ 
			eff_file >> values[i];
		}

		energy.push_back(values[0]);
		efficiency.push_back(values[1]);
		if(eff_file.eof()){ break; }
	}	
	eff_file.close();

	if(energy.size() != efficiency.size()){ return false; }
	return true;
}

Efficiency::Efficiency(){
	small_energy = NULL; small_efficiency = NULL;
	med_energy = NULL; med_efficiency = NULL;
	large_energy = NULL; large_efficiency = NULL;
	NsmallEff = 0; NmedEff = 0; NlargeEff = 0;
	init_small = false; 
	init_med = false; 
	init_large = false;
}

Efficiency::~Efficiency(){
	if(init_small){ delete[] small_energy; delete[] small_efficiency; }
	if(init_med){ delete[] med_energy; delete[] med_efficiency; }
	if(init_large){ delete[] large_energy; delete[] large_efficiency; }
}

// Load small bar efficiency data
unsigned int Efficiency::ReadSmall(const char* fname){
	std::vector<double> E, Eff;
	if(init_small || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	small_energy = new double[E.size()];
	small_efficiency = new double[E.size()];
	for(unsigned int i = 0; i < E.size(); i++){
		small_energy[i] = E[i];
		small_efficiency[i] = Eff[i];
	}
	init_small = true;

	NsmallEff = E.size();
	return NsmallEff;
}

// Load medium bar efficiency data
unsigned int Efficiency::ReadMedium(const char* fname){
	std::vector<double> E, Eff;
	if(init_med || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	med_energy = new double[E.size()];
	med_efficiency = new double[E.size()];
	for(unsigned int i = 0; i < E.size(); i++){
		med_energy[i] = E[i];
		med_efficiency[i] = Eff[i];
	}
	init_med = true;

	NmedEff = E.size();
	return NmedEff;
}

// Load medium bar efficiency data
unsigned int Efficiency::ReadLarge(const char* fname){
	std::vector<double> E, Eff;
	if(init_large || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	large_energy = new double[E.size()];
	large_efficiency = new double[E.size()];
	for(unsigned int i = 0; i < E.size(); i++){
		large_energy[i] = E[i];
		large_efficiency[i] = Eff[i];
	}
	init_large = true;

	NlargeEff = E.size();
	return NlargeEff;
}

// Return the interpolated value for the input energy
double Efficiency::GetSmallEfficiency(double Energy){
	if(!init_small || NsmallEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < small_energy[0]){ efficiency = small_efficiency[0]; }
	else if(Energy > small_energy[NsmallEff-1]){ efficiency = small_efficiency[NsmallEff-1]; }
	else{
		for(unsigned int i = 1; i < NsmallEff; i++){
			if(Energy >= small_energy[i-1] && Energy <= small_energy[i]){
				efficiency = Interpolate(small_energy[i-1],small_efficiency[i-1],
							 small_energy[i],small_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
}

// Return the interpolated value for the input energy
double Efficiency::GetMediumEfficiency(double Energy){
	if(!init_med || NmedEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < med_energy[0]){ efficiency = med_efficiency[0]; }
	else if(Energy > med_energy[NmedEff-1]){ efficiency = med_efficiency[NmedEff-1]; }
	else{
		for(unsigned int i = 1; i < NmedEff; i++){
			if(Energy >= med_energy[i-1] && Energy <= med_energy[i]){
				efficiency = Interpolate(med_energy[i-1],med_efficiency[i-1],
							 med_energy[i],med_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
}	

// Return the interpolated value for the input energy
double Efficiency::GetLargeEfficiency(double Energy){
	if(!init_large || NlargeEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < large_energy[0]){ efficiency = large_efficiency[0]; }
	else if(Energy > large_energy[NlargeEff-1]){ efficiency = large_efficiency[NlargeEff-1]; }
	else{
		for(unsigned int i = 1; i < NlargeEff; i++){
			if(Energy >= large_energy[i-1] && Energy <= large_energy[i]){
				efficiency = Interpolate(large_energy[i-1],large_efficiency[i-1],
							 large_energy[i],large_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
}

/////////////////////////////////////////////////////////////////////
// Material
/////////////////////////////////////////////////////////////////////

void Material::_initialize(){
	vikar_name = "unknown";
	density = 1.0;
	rad_length = 0.0;
	avgZ = 0.0;
	avgA = 0.0;
	lnIbar = 0.0;
	num_elements = 0;
	total_elements = 0;
	num_per_molecule = NULL;
	element_Z = NULL;
	element_A = NULL;
	element_I = NULL;
	init = false;
}

// Calculate the average atomic charge, mass, and ionization potential for the material
void Material::_calculate(){
	double numerator1 = 0.0;
	double numerator2 = 0.0;
	double numerator3 = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		numerator1 += num_per_molecule[i]*element_Z[i]; // Average Z
		numerator2 += num_per_molecule[i]*element_A[i]; // Average A
		numerator3 += num_per_molecule[i]*element_Z[i]*element_I[i]; // Average I
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		denominator += num_per_molecule[i]; // Total number of elements
	}
	avgZ = numerator1/denominator;
	avgA = numerator2/denominator;
	lnIbar = numerator3/numerator1;
	edens = avagadro*avgZ*(density*(1E6))/avgA;

	// Calculate the radiation length for the target material
	// for use in angular straggling calculations later.
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
	rad_length = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){ 
		rad_length += ((element_A[i]*num_per_molecule[i]/avgA) / radlength(element_A[i],element_Z[i])); 
	}

	rad_length = 1.0/rad_length; 
}

Material::Material(){
	_initialize();
}

Material::Material(unsigned int num_elements_){
	_initialize();
	Init(num_elements_);
}

Material::~Material(){
	if(init){
		delete[] num_per_molecule;
		delete[] element_Z;
		delete[] element_A;
		delete[] element_I;
	}
}

bool Material::Init(unsigned int num_elements_){
	if(init){ return false; }
	num_elements = num_elements_;
	num_per_molecule = new unsigned int[num_elements];
	element_Z = new double[num_elements];
	element_A = new double[num_elements];
	element_I = new double[num_elements];
	for(unsigned int i = 0; i < num_elements; i++){
		num_per_molecule[i] = 0;
		element_Z[i] = 0.0;
		element_A[i] = 0.0;
		element_I[i] = 0.0;
	}

	rad_length = 0.0;
		
	init = true;
	return true;
}

void Material::SetElements(unsigned int *num_per_molecule_, double *Z_, double *A_){
	for(unsigned int i = 0; i < num_elements; i++){
		num_per_molecule[i] = num_per_molecule_[i];
		total_elements += num_per_molecule_[i];
		element_Z[i] = Z_[i];
		element_A[i] = A_[i];
		element_I[i] = GetIonPot(element_Z[i]);
	}
	
	_calculate();
}

bool Material::ReadMatFile(const char* filename_){
	if(init){ return false; }
	std::ifstream input_file(filename_);
	if(!input_file.good()){ return false; }

	unsigned int count = 0;
	std::string line;
	while(true){
		getline(input_file, line);
		if(input_file.eof()){ break; }
		if(line[0] == '#'){ continue; } // Commented line
		line = Parse(line);
		
		if(count == 0){ vikar_name = line; }
		else if(count == 1){ density = (double)atof(line.c_str()); }
		else if(count == 2){ 
			num_elements = (unsigned int)atoi(line.c_str());
			Init(num_elements);
			
			for(unsigned int i = 0; i < num_elements; i++){
				getline(input_file, line); line = Parse(line); element_Z[i] = (double)atof(line.c_str()); 
				getline(input_file, line); line = Parse(line); element_A[i] = (double)atof(line.c_str());  
				getline(input_file, line); line = Parse(line); num_per_molecule[i] = (unsigned int)atoi(line.c_str()); 
				element_I[i] = GetIonPot(element_Z[i]);
			}
			
			_calculate();
		}
		count++;
	}
	if(count <= 2){ return false; }

	return true;
}

// Return the ionization potential (MeV) for a particle with a given atomic charge Z_
// Z_ is the atomic number of the ion of interest
double Material::GetIonPot(unsigned int Z_){ 
	if(Z_ <= 13){ return ionpot[Z_-1]*(1E-6); }
	return((9.76*Z_ + 58.8*pow(1.0*Z_, -0.19))*(1E-6));
}

/*double shell(double e){
	//          to calculate the shell correction term in stopping power
	//     calculations(spar-armstrong & chandler-ornl-4869 [1973])
	//
	//     input - e energy in mev
	//             avip mean ionization potential of stopping medium (mev)
	//             avz mean atomic number of stopping medium
	 
	double be2; 
	double gnu2; 
	double f1, f2, xl, xl1, xl2; 
	double output, x, xlog;
	
	const double a1 = 0.422377E-6, a2 = 3.858019E-9, b1 = 0.304043E-7; 
	const double b2 = -0.1667989E-9, c1 = -0.38106E-9, c2 = 0.157955E-11; 
	const double emp = 938.25920, zal = 12.950, zh2o = 3.340; 
	const double emp1 = 0.00106580356; 
	
	const double p1 = 4.774248E-4, p2 = 1.143478E-3, p3 = -5.63392E-2; 
	const double p4 = 4.763953E-1, p5 = 4.844536E-1; 
	const double w1 = -1.819954E-6, w2 = -2.232760E-5, w3 = 1.219912E-4; 
	const double w4 = 1.837873E-3, w5 = -4.457574E-3, w6 = -6.837103E-2; 
	const double w7 = 5.266586E-1, w8 = 3.743715E-1; 
	const double const1 = 0.021769515; 

	if(e >= 8.0){
		gnu2 = 1.0/((e*emp1)*(e*emp1+2.00)); 
		f1 =  gnu2*(a1+gnu2*(b1+c1*gnu2)); 
		f2 =  gnu2*(a2+gnu2*(b2+c2*gnu2)); 
		return avip*avip*(f1 + f2*avip)/avz; 
	} 
	
	be2 = std::sqrt(2*e/emp); 
	output = const1 + std::log(be2) - std::log(avip); 
	x = 18769.00*be2/avz; 
	xlog=std::log(x);
	xl1 = 0.0;
	
	if(avz > zal){
		xl = p5 + xlog*(p4+xlog*(p3+xlog*(p2+xlog*p1)));
		xl = std::exp(xl);
		if(avz > zal){ output = output - xl; }
		else{
			xl2 = xl;
			xl = xl1 + (avz-zh2o)/(zal-zh2o)*(xl2-xl1);
			output = output - xl;
		}
	}
	else{
		xl1 = w8 + xlog*(w7+xlog*(w6+xlog*(w5+xlog*(w4+xlog*(w3+xlog*(w2+xlog*w1))))));
		xl1 = std::exp(xl1);
		if(avz > zh2o){
			xl = p5 + xlog*(p4+xlog*(p3+xlog*(p2+xlog*p1)));
			xl = std::exp(xl);
			if(avz > zal){ output = output - xl; }
			else{
				xl2 = xl;
				xl = xl1 + (avz-zh2o)/(zal-zh2o)*(xl2-xl1);
				output = output - xl;
			}	
		}
		else{
			xl = xl1;
			output = output - xl;
		}
	} 
	
	return output;
} */

// Return the shell correction term for a particle with a given energy_ in this material
// energy_ in MeV
double Material::ShellCorrect(double energy_){
	if(energy_ >= 8.0){
		double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
		double f1 = ((a1/nu2) + (b1/(nu2*nu2)) + (c1*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		double f2 = ((a2/nu2) + (b2/(nu2*nu2)) + (c2*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		double Ibar = std::exp(lnIbar); 
		return ((f1*pow(Ibar, 2.0) + f2*pow(Ibar, 3.0))/GetAverageZ());
	}
	else{
		// Aluminum and water stopping power
		double beta2 = Beta2(energy_, 1.0);
		double Zbar = GetAverageZ();
		double x = 18769.0*beta2/Zbar; // 137^2 = 187690
		double output = std::log(2*electron_RME*beta2) - lnIbar;
		
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
double Material::DensityEffect(double energy_){
	double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
	double output = std::log(alpha) + std::log(edens) + 2.0*(std::log(nu2) - lnIbar) - 1;
	if(output < 0.0){ return 0.0; }
	return output;
}

// Return the stopping power for a proton with a given energy_ in this material
// energy_ in MeV
double Material::Pstop(double energy_){
	double beta2 = Beta2(energy_, proton_RME);
	double beta = std::sqrt(beta2);
	double output = edens*coeff*pow(Zeff(beta, 1.0), 2.0)/beta2; // C^4/(MeV*cm^3)
	output *= -1.0*(std::log(2.0*electron_RME*(beta2/(1-beta2)))-beta2-lnIbar-ShellCorrect(energy_)-0.5*DensityEffect(energy_)); // 
	return output;
}

// Return the stopping power (MeV/m) for a particle with given energy_, Z_, and mass_ in this material
// energy_ in MeV and mass_ in MeV/c^2
double Material::StopPower(double energy_, double Z_, double mass_){
	double beta = std::sqrt(Beta2(energy_, mass_));
	double output = pow((Zeff(beta, Z_)/Zeff(beta, 1.0)), 2.0);
	output *= Pstop(energy_*(proton_RME/mass_));
	return output;
}

// Return the range (m) for a particle with given energy_, Z_, and mass_, in this material
// energy_ in MeV and mass_ in MeV/c^2
double Material::Range(double energy_, double Z_, double mass_){
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

void Material::Print(){
	std::cout << " Name: " << vikar_name << std::endl;
	std::cout << "  Number of unique elements: " << num_elements << std::endl;
	for(unsigned int i = 0; i < num_elements; i++){
		std::cout << "   Element " << i+1 << ": N = " << num_per_molecule[i] << ", Z = " << element_Z[i];
		std::cout << ", A = " << element_A[i] << ", I = " << element_I[i] << std::endl;
	}
	std::cout << "  Average Z: " << avgZ << "\n";
	std::cout << "  Average A: " << avgA << "\n";
	std::cout << "  Density: " << density << " g/cm^3\n";
	std::cout << "  Edensity: " << edens << " 1/m^3\n";
	std::cout << "  Rad Length: " << rad_length << " mg/cm^2\n";
	std::cout << "  lnIbar: " << lnIbar << "\n";
}

/////////////////////////////////////////////////////////////////////
// Target
/////////////////////////////////////////////////////////////////////

Target::Target(){
	thickness = 0.0;
	Zthickness = 0.0;
	density = 1.0;
	rad_length = 0.0;
	angle = 0.0;
	physical = new Planar();
}

Target::Target(unsigned int num_elements_){
	thickness = 0.0;
	Zthickness = 0.0;
	density = 1.0;
	rad_length = 0.0;
	angle = 0.0;
	physical = new Planar();
}

void Target::SetThickness(double thickness_){ 
	thickness = thickness_;	
	physical->SetSize(1.0, 1.0, (thickness/density)*1E-5); // Only the thickness matters
}

void Target::SetAngle(double angle_){
	angle = angle_;	
	Zthickness = dabs(thickness / std::cos(angle));
	physical->SetRotation(angle, 0.0, 0.0);
}

// Get the depth into the target at which the reaction occurs
// offset_ is the global position where the beam particle originates
// direction_ is the direction of the beam particle entering the target
// intersect is the global position where the beam particle intersects the front face of the target
// interact is the global position where the beam particle reacts inside the target
double Target::GetInteractionDepth(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, Vector3 &interact){
	Vector3 temp;
	double zdist = physical->GetApparentThickness(offset_, direction_, 0, 2, intersect, temp); // The target thickness the ray sees
	if(thickness == -1){ 
		std::cout << " Beam does not travel through target!\n"; 
		return -1;
	}
	
	zdist *= frand(); // Reaction occurs at a random depth into the target
	interact = intersect + direction_*zdist;
	return zdist; 
}

// Determine the new direction of a particle inside the target due to angular straggling
// direction_ and new_direction have x,y,z format and are measured in meters
bool Target::AngleStraggling(const Vector3 &direction_, double A_, double Z_, double E_, Vector3 &new_direction){
	// strag_targ 1.0 written by S.D.Pain on 20/01/2004
	//
	// strag_targ 1.1 modified by S.D.Pain on 7/03/2005
	// to untilise transform2.0
	//
	// Subroutine to calculate the angular straggling of an ion in
	// the target. The A,Z of the ion are read in, along with the
	// theta,phi and energy of the ion. The average radiation length
	// of the target material is read in as X.
	//
	// The calculation of the width of the scattering distribution
	// is calculated by straggleA. A Gaussian weighted scattering angle
	// based on this width is calculated, using rndgauss0
	//
	// The new theta,phi to which the ion is scattered is returned.
		
	double theta_scatW; 
	double theta_scat, phi_scat; 
	
	// Calculate the straggling width
	straggleA(theta_scatW, E_, Z_, A_, thickness, rad_length); 
	
	// Select the scattering angle of the ion wrt its initial direction
	theta_scat = rndgauss0(theta_scatW); 
	theta_scat = std::sqrt(pow(theta_scat, 2.0)*2.0); 
	phi_scat = frand()*2.0*pi; 
	
	// Determine the absolute lab angle to which the ion is scattered
	Matrix3 matrix(theta_scat, phi_scat);
	new_direction = direction_;
	matrix.Transform(new_direction);
	return true;
}
