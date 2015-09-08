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
const double e = 1.602176565e-19; // C
const double h = 4.135667516e-21; // MeV * s
const double coeff = 4*pi*pow(bohr_e_radius, 2.0)*electron_RME; // m^2 * MeV (checked)
const double alpha = h*h*9E16*bohr_e_radius/pi; // MeV^2 * m^3
const double avagadro = 6.0221413e+23; // 1/mol
const double amu2mev = 931.494061; // MeV

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

bool RangeTable::Init(unsigned int num_entries_, double startE_, double stopE_, double Z_, double mass_, Material *mat_){
	if(!_initialize(num_entries_)){ return false; }
	
	// Use Material to fill the arrays
	double step = (stopE_-startE_)/(num_entries_-1);
	for(unsigned int i = 0; i < num_entries_; i++){
		energy[i] = startE_ + i*step; // Energy in MeV
		range[i] = mat_->Range(energy[i], Z_, mass_); // Range in m
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

double RangeTable::GetNewE(double energy_, double dist_){
	if(!use_table){ return -1; }
	double r1 = GetRange(energy_);
	if(r1 - dist_ > 0.0){ return GetEnergy(r1 - dist_); } // The particle loses some energy in the material
	else{ return 0.0; } // The particle stops in the material
	return -1;
}

void RangeTable::Print(){
	if(!use_table){ return; }
	for(unsigned int i = 0; i < num_entries; i++){
		std::cout << i << "\t" << energy[i] << "\t" << range[i] << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////
// Efficiency
/////////////////////////////////////////////////////////////////////

// Read NewVIKAR VANDLE efficiency file and load data into arrays
// Returns the number of data points which are found in the file
// Assumes the following efficiency file format for each point in file
// Lab__energy(MeV) Efficiency
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
	Mmass = 0.0;
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
	use_eloss = true;
}

// Calculate the average atomic charge, mass, and ionization potential for the material
void Material::_calculate(){
	double total_num_dens = avagadro*(density*(1E6))/Mmass; // The number density of the material (1/m^3)
	double num_dens; // The number density of each element in the material (1/m^3)
	double ztot = 0.0; // Z total
	double atot = 0.0; // A total
	lnIbar = 0.0;
	edens = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		ztot += num_per_molecule[i]*element_Z[i];
		atot += num_per_molecule[i]*element_A[i];
		num_dens = num_per_molecule[i]*total_num_dens;
		
		lnIbar += num_dens*element_Z[i]*element_I[i];
		edens += num_dens*element_Z[i];
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		denominator += num_per_molecule[i]; // Total number of elements
	}
	avgZ = ztot/denominator;
	avgA = atot/denominator;
	lnIbar = lnIbar/edens;

	// Calculate the radiation length for the target material
	// for use in angular straggling calculations later.
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
	rad_length = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){ 
		rad_length += ((element_A[i]*num_per_molecule[i]/avgA) / radlength(element_A[i],element_Z[i])); 
	}

	rad_length = 1.0/rad_length; 
}

// Return the ionization potential (MeV) for a particle with a given atomic charge Z_
// Z_ is the atomic number of the ion of interest
double Material::_ionpot(unsigned int Z_){ 
	if(Z_ <= 13){ return ionpot[Z_-1]*(1E-6); }
	return((9.76*Z_ + 58.8*pow(1.0*Z_, -0.19))*(1E-6));
}

// Return the shell correction term for a proton with a given KE (energy_) in this material
// This corrects the assumption that the proton's velocity is much greater than the bound electron velocity
// energy_ in MeV
double Material::_shell(double energy_){
	double Ibar = std::exp(lnIbar);
	if(energy_ >= 8.0){
		double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
		double f1 = ((a1/nu2) + (b1/(nu2*nu2)) + (c1*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		double f2 = ((a2/nu2) + (b2/(nu2*nu2)) + (c2*std::sqrt(nu2)/(nu2*nu2*nu2)))*(1E-6); // MeV
		return ((f1*pow(Ibar, 2.0) + f2*pow(Ibar, 3.0))/avgZ);
	}
	else{
		return 0.0;
		// Aluminum and water stopping power
		/*double beta2 = _b2(energy_, proton_RME);
		double x = 18769.0*beta2/avgZ; // 137^2 = 187690
		double lnx = std::log(x);
		double output = std::log(2*electron_RME*beta2) - lnIbar;
		double xl, xl1;
	
		const double p1 = 4.774248E-4, p2 = 1.143478E-3, p3 = -5.63392E-2; 
		const double p4 = 4.763953E-1, p5 = 4.844536E-1; 
		const double w1 = -1.819954E-6, w2 = -2.232760E-5, w3 = 1.219912E-4; 
		const double w4 = 1.837873E-3, w5 = -4.457574E-3, w6 = -6.837103E-2; 
		const double w7 = 5.266586E-1, w8 = 3.743715E-1; 
		
		// Need to calculate L(x)
		if(avgZ >= 13.0){ // Use L2(x) (proton stopping in Al)
			xl = std::exp(p5 + lnx*(p4+lnx*(p3+lnx*(p2+lnx*p1))));
			if(avgZ > 13.0){ output = output - xl; }
			else{ output = output - (avgZ-3.33)/(13.0-3.33)*xl; }
		}
		else{
			xl1 = std::exp(w8 + lnx*(w7+lnx*(w6+lnx*(w5+lnx*(w4+lnx*(w3+lnx*(w2+lnx*w1)))))));
			if(avgZ <= 3.33){ // Use L1(x) (proton stopping in H20)
				xl = std::exp(p5 + lnx*(p4+lnx*(p3+lnx*(p2+lnx*p1))));
				if(avgZ > 13.0){ output = output - xl; }
				else{ output = output - xl1 + (avgZ-3.33)/(13.0-3.33)*(xl-xl1); }
			}
			else{ // Interpolate between L1 and L2
				output = output - std::exp(w8 + lnx*(w7+lnx*(w6+lnx*(w5+lnx*(w4+lnx*(w3+lnx*(w2+lnx*w1)))))));
			}
		}
		return output;*/
	}
	return 0.0;
}

// Return the density effect correction term (unitless) for a proton with a given KE (energy_) in this material
// This corrects for dielectric effects in the medium
// energy_ in MeV
double Material::_density(double energy_){
	double nu2 = (energy_/proton_RME)*((energy_/proton_RME) + 2); // unitless
	double output = std::log(alpha*edens*nu2/(lnIbar*lnIbar)) - 1;
	if(output < 0.0){ return 0.0; }
	return output;
}

// Return the stopping power for a proton with a given KE (energy_) in this material
// energy_ in MeV
double Material::_pstop(double energy_){
	double beta2 = _b2(energy_, proton_RME);
	double beta = std::sqrt(beta2);
	double output = edens*coeff*pow(_zeff(beta, 1.0), 2.0)/beta2; // Stopping coefficient
	double L0 = std::log(2.0*electron_RME*(beta2/(1-beta2)))-beta2-lnIbar; // Primary stopping number
	//double L0 = std::log(2.0*electron_RME*(beta2/(1-beta2)))-beta2-lnIbar-ShellCorrect(energy_)-0.5*DensityEffect(energy_); // Primary stopping number
	return -2.5*output*(L0);
}

// Return the range (m) for a proton with given KE (energy_) in this material
// energy_ in MeV
double Material::_prange(double energy_, unsigned int num_iterations_){
	double sum = 0.0;
	double e1, e2;
	double step = energy_/num_iterations_;
	for(unsigned int i = 1; i < num_iterations_; i++){
		e1 = step*i; e2 = step*(i+1);
		sum += 0.5*(1.0/_pstop(e2)+1.0/_pstop(e1))*step;
	}

	return sum;
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
		element_I[i] = _ionpot(element_Z[i]);
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
		else if(count == 2){ Mmass = (double)atof(line.c_str()); }
		else if(count == 3){  
			num_elements = (unsigned int)atoi(line.c_str());
			Init(num_elements);
			
			for(unsigned int i = 0; i < num_elements; i++){
				getline(input_file, line); line = Parse(line); element_Z[i] = (double)atof(line.c_str()); 
				getline(input_file, line); line = Parse(line); element_A[i] = (double)atof(line.c_str());  
				getline(input_file, line); line = Parse(line); num_per_molecule[i] = (unsigned int)atoi(line.c_str()); 
				element_I[i] = _ionpot(element_Z[i]);
			}
			
			_calculate();
		}
		count++;
	}
	if(count <= 3){ return false; }

	return true;
}

// Return the stopping power (MeV/m) for a particle with given KE (energy_), Z_, and mass_ in this material
// energy_ in MeV and mass_ in MeV/c^2
double Material::StopPower(double energy_, double Z_, double mass_){
	double beta = std::sqrt(_b2(energy_, mass_));
	double output = pow((_zeff(beta, Z_)/_zeff(beta, 1.0)), 2.0);
	output *= _pstop(energy_*(proton_RME/mass_));
	return output;
}

// Return the range (m) for a particle with given KE (energy_), Z_, and mass_, in this material
// energy_ in MeV and mass_ in MeV/c^2
double Material::Range(double energy_, double Z_, double mass_){
	unsigned int iterations = (unsigned int)Order(energy_)*10;
	double Eci = _energy(0.04*pow(Z_, 2.0/3.0), mass_);
	double sum = 0.0;
	if(energy_ <= Eci){
		double e1, e2;
		double step = energy_/iterations;
		for(unsigned int i = 1; i < iterations; i++){
			e1 = step*i; e2 = step*(i+1);
			sum += 0.5*(1.0/StopPower(e2, Z_, mass_)+1.0/StopPower(e1, Z_, mass_))*step;
		}
	}
	else{
		double e1, e2;
		double step = Eci/iterations;
		
		// Ri(Eci, Zi, mi)
		for(unsigned int i = 1; i < iterations; i++){
			e1 = step*i; e2 = step*(i+1);
			sum += 0.5*(1.0/StopPower(e2, Z_, mass_)+1.0/StopPower(e1, Z_, mass_))*step;
		}
		
		sum += (mass_/(proton_RME*Z_*Z_))*(_prange(energy_*proton_RME/mass_, iterations) - _prange(Eci*proton_RME/mass_, iterations));
	}
	return sum;
}

	// Use Birks' equation to calculate the light output for a particle at a given energy_ in this material
	// energy_ in MeV, mass_ in MeV/c^2, L0_ in 1/MeV, kB_ in m/MeV, and C_ in (m/MeV)^2
double Material::Birks(double energy_, double Z_, double mass_, double L0_, double kB_, double C_){
	unsigned int iterations = (unsigned int)Order(energy_)*10;
	double sum = 0.0;
	double dedx1, dedx2;
	double step = energy_/iterations;
	for(unsigned int i = 1; i < iterations; i++){
		dedx1 = StopPower(step*i, Z_, mass_);
		dedx2 = StopPower(step*(i+1), Z_, mass_);
		sum += 0.5*L0_*((1.0/(1.0 + kB_*dedx1 + C_*dedx1*dedx1)) + (1.0/(1.0 + kB_*dedx2 + C_*dedx2*dedx2)))*step;
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
	std::cout << "  Molar Mass: " << Mmass << " g/mol\n";
	std::cout << "  Rad Length: " << rad_length << " mg/cm^2\n";
	std::cout << "  lnIbar: " << lnIbar << "\n";
}

void Material::Print(std::ofstream *file_){
	if(!file_ || !file_->good()){ return; }
	(*file_) << "Name\t" << vikar_name << std::endl;
	(*file_) << "NumElements\t" << num_elements << std::endl;
	for(unsigned int i = 0; i < num_elements; i++){
		(*file_) << "Element" << i+1 << "\tN=" << num_per_molecule[i] << ", Z=" << element_Z[i];
		(*file_) << ", A=" << element_A[i] << ", I=" << element_I[i] << std::endl;
	}
	(*file_) << "AvgZ\t" << avgZ << "\n";
	(*file_) << "AvgA\t" << avgA << "\n";
	(*file_) << "Density\t" << density << " g/cm^3\n";
	(*file_) << "eDensity\t" << edens << " 1/m^3\n";
	(*file_) << "MolarMass\t" << Mmass << " g/mol\n";
	(*file_) << "RadLength\t" << rad_length << " mg/cm^2\n";
	(*file_) << "lnIbar\t" << lnIbar << "\n";
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
	physical->SetSize(1.0, 1.0, thickness/(density*1E5)); // Only the thickness matters
}

void Target::SetAngle(double angle_){
	angle = angle_;	
	Zthickness = dabs(thickness / std::cos(angle));
	physical->SetRotation(angle, 0.0, 0.0);
}

void Target::SetDensity(double density_){
	density = density_;
	
	// The density changes the physical thickness of the detector
	// so we need to update the thickness
	physical->SetSize(1.0, 1.0, thickness/(density*1E5));
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

/////////////////////////////////////////////////////////////////////
// Particle
/////////////////////////////////////////////////////////////////////

bool Particle::SetMaterial(Material *mat_, double Ebeam_, double Espread_){
	if(Z <= 0.0 || init){ return false; }
	mat = mat_;
	maxE = Ebeam_ + 2*Espread_;
	
	if(table.Init(100, 0.1, maxE, Z, mass, mat)){
		init = true;
		return true;
	}
	return false;
}
