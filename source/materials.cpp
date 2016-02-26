/** \file materials.cpp
 * \brief Handles particles, isotopes, and their interactions w/ matter.
 *
 * This file contains classes to handle particle/isotope interaction with
 * matter, detector efficiencies, particle ranges, etc.
 *
 * \author C. R. Thornsberry
 * \date Feb. 26th, 2016
 */
#include "vandmc_core.h"
#include "materials.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals 
/////////////////////////////////////////////////////////////////////

// Constants for use in stopping power calcuations
const double electron_RME = 0.510998928; // MeV
const double proton_RME = 938.272046; // MeV
const double neutron_RME = 939.565378; // MeV
const double bohr_e_radius = 2.817940326836615e-15; // m
const double e_charge = 1.60217662E-19; // C
const double avagadro = 6.0221413e+23; // 1/mol
const double mev2amu = 1.0/931.494061; // (amu*c^2)/MeV
const double mev2kg = 1.783E-30; // (kg*c^2)/MeV

// Ionization potentials for Z=1 to Z=100 (in eV).
const float potentials[100] = {19.2, 41.8, 40, 63.7, 76, 78, 82, 95, 115, 137, 
	                           149, 156, 166, 173, 173, 180, 174, 188, 190, 191, 
	                           216, 233, 245, 257, 272, 286, 297, 311, 322, 330, 
	                           334, 350, 347, 348, 357, 352, 363, 366, 379, 393, 
	                           417, 424, 428, 441, 449, 470, 470, 469, 488, 488, 
	                           487, 485, 491, 482, 488, 491, 501, 523, 535, 546, 
	                           560, 574, 580, 591, 614, 628, 650, 658, 674, 684, 
	                           694, 705, 718, 727, 736, 746, 757, 790, 790, 800, 
	                           810, 823, 823, 830, 825, 794, 827, 826, 841, 847, 
	                           878, 890, 902, 921, 934, 939, 952, 966, 980, 994};

/////////////////////////////////////////////////////////////////////
// RangeTable
/////////////////////////////////////////////////////////////////////

/// Initialize range table arrays.
bool RangeTable::_initialize(const unsigned int &num_entries_){
	if(use_table){ return false; }
	num_entries = num_entries_;
	energy = new double[num_entries];
	dedx = new double[num_entries];
	range = new double[num_entries];
	use_table = true;
	return true;
}

/// Interpolate between two points
double RangeTable::_interpolate(double *x_, double *y_, const double &val_){
	if(!x_ || !y_){ return -1; }
	else if(val_ < x_[0]){ return 0.0; }
	for(unsigned int i = 0; i < num_entries-1; i++){
		if(val_ == x_[i]){ return y_[i]; }
		else if(val_ == x_[i+1]){ return y_[i+1]; }
		else if(val_ >= x_[i] && val_ <= x_[i+1]){
			// Interpolate and return the result
			return (((y_[i+1]-y_[i])/(x_[i+1]-x_[i]))*(val_-x_[i])+y_[i]);
		}
	}
	return -1;
}

RangeTable::RangeTable(){ 
	use_table = false; 
	use_birks = false;
}

/// Constructor to set the number of table entries.
RangeTable::RangeTable(const unsigned int &num_entries_){
	_initialize(num_entries_);
	use_birks = false;
}

/// Destructor.
RangeTable::~RangeTable(){
	if(use_table){
		delete[] energy;
		delete[] dedx;
		delete[] range;
	}
	if(use_birks){ delete[] birks; }
}
	
/// Initialize arrays but do not fill them.
bool RangeTable::Init(const unsigned int &num_entries_){
	return _initialize(num_entries_);
}

/// Initialize arrays and fill them using Material.
bool RangeTable::Init(const unsigned int &num_entries_, const double &startE_, const double &stopE_, const double &Z_, const double &mass_, Material *mat_){
	if(!_initialize(num_entries_)){ return false; }
	
	// Use Material to fill the stopping power arrays.
	step = (stopE_-startE_)/(num_entries_-1);
	for(unsigned int i = 0; i < num_entries_; i++){
		energy[i] = startE_ + i*step; // Energy in MeV
		dedx[i] = -1*mat_->StopPower(energy[i], Z_, mass_); // Stopping power in MeV/m;
	}
	
	range[0] = 0.0;
	
	// Calculate ranges.
	for(unsigned int i = 1; i < num_entries_; i++){
		range[i] = range[i-1] - 0.5*(1.0/dedx[i-1] + 1.0/dedx[i])*step;
	}
	
	return true;
}

/// Initialize the birks light response array.
bool RangeTable::InitBirks(double L0_, double kB_, double C_/*=0.0*/){
	if(!use_table || use_birks){ return false; }
	
	birks = new double[num_entries];
	birks[0] = 0.0;
	
	for(unsigned int i = 1; i < num_entries; i++){
		birks[i] = birks[i-1] + L0_*0.5*(1.0/(1.0 + kB_*dedx[i-1] + C_*dedx[i-1]*dedx[i-1]) + 1.0/(1.0 + kB_*dedx[i] + C_*dedx[i]*dedx[i]))*step;
	}
	
	use_birks = true;
	
	return true;
}

/// Manually set a data point with an energy and a range.
bool RangeTable::Set(const unsigned int &pt_, const double &energy_, const double &range_){
	if(!use_table){ return false; }
	if(pt_ > num_entries){ return false; }
	energy[pt_] = energy_;
	range[pt_] = range_;
	return true;
}

/// Get the particle range at a given energy using linear interpolation.
double RangeTable::GetRange(const double &energy_){
	if(!use_table){ return -1; }
	return _interpolate(this->energy, this->range, energy_);
}

/// Get the particle energy at a given range using linear interpolation.
double RangeTable::GetEnergy(const double &range_){
	if(!use_table){ return -1; }
	return _interpolate(this->range, this->energy, range_);
}

/// Get the scintillator light response due to a particle traversing a material with given kinetic energy.
double RangeTable::GetLRfromKE(const double &energy_){
	if(!use_birks){ return -1; }
	return _interpolate(this->energy, this->birks, energy_);
}

/// Get the kinetic energy of a particle which produces a given light response in a scintillator.
double RangeTable::GetKEfromLR(const double &response_){
	if(!use_birks){ return -1; }
	return _interpolate(this->birks, this->energy, response_);
}

/// Get the new energy of a particle traversing a distance through a material.
double RangeTable::GetNewE(const double &energy_, const double &dist_){
	double dummy;
	return GetNewE(energy_, dist_, dummy);
}

/// Get the new energy of a particle traversing a distance through a material.
double RangeTable::GetNewE(const double &energy_, const double &dist_, double &dist_traveled){
	if(!use_table){ return -1; }
	dist_traveled = GetRange(energy_);
	if(dist_traveled < 0.0){ return -1; }
	if(dist_traveled - dist_ > 0.0){ // The particle loses some energy in the material
		dist_traveled = dist_;
		return GetEnergy(dist_traveled - dist_); 
	} 
	else{ return 0.0; } // The particle stops in the material
	return -1;
}

/// Return the range and energy for an entry in the table.
bool RangeTable::GetEntry(const unsigned int &entry_, double &E, double &R){
	if(!use_table || entry_ >= num_entries){ return false; }
	E = energy[entry_]; R = range[entry_];
	return true;
}

/// Print range table entries to the screen.
void RangeTable::Print(){
	if(!use_table){ return; }
	for(unsigned int i = 0; i < num_entries; i++){
		std::cout << i << "\t" << energy[i] << "\t" << dedx[i] << "\t" << range[i];
		if(use_birks){ std::cout << "\t" << birks[i]; }
		std::cout << std::endl;
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

/** Initialize all variables with default values.
  */
void Material::_initialize(){
	vikar_name = "unknown";
	density = 1.0;
	Mmass = 0.0;
	rad_length = 0.0;
	avgZ = 0.0;
	avgA = 0.0;
	lnIbar = 0.0;
	coeff = 0.0;
	num_elements = 0;
	total_elements = 0;
	num_per_molecule = NULL;
	element_Z = NULL;
	element_A = NULL;
	element_I = NULL;
	weight = NULL;
	init = false;
	use_eloss = true;
}

/** Calculate the average atomic charge, mass, and ionization potential for the material.
  */
void Material::_calculate(){
	// Calculate the total molar mass of the material (g/mol).
	Mmass = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		Mmass += num_per_molecule[i] * element_A[i];
	}
	
	// Calculate the fractional weight for each element.
	for(unsigned int i = 0; i < num_elements; i++){
		weight[i] = num_per_molecule[i] * element_A[i] / Mmass;
	}

	double ztot = 0.0; // Z total
	double atot = 0.0; // A total
	lnIbar = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		ztot += num_per_molecule[i]*element_Z[i];
		atot += num_per_molecule[i]*element_A[i];
		
		lnIbar += weight[i]*(element_Z[i]/element_A[i])*std::log(element_I[i]);
	}
	double denominator = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		denominator += num_per_molecule[i]; // Total number of elements
	}
	avgZ = ztot/denominator;
	avgA = atot/denominator;

	denominator = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		denominator += weight[i] * element_Z[i] / element_A[i];
	}
	
	lnIbar = lnIbar/denominator;

	coeff = 4.0 * pi * avagadro * std::pow(bohr_e_radius, 2.0) * electron_RME * (density * 1E6); // (MeV * g / (mol * m))

	// Calculate the radiation length for the target material
	// for use in angular straggling calculations later.
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
	rad_length = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){ 
		rad_length += ((element_A[i]*num_per_molecule[i]/avgA) / radlength(element_A[i],element_Z[i])); 
	}

	rad_length = 1.0/rad_length; 
}

/** Return the velocity of a particle with a given energy and mass relative to c.
  * \param[in] energy_ The kinetic energy of the incoming particle (MeV).
  * \param[in] mass_ The mass of the incoming particle (MeV/c^2).
  * \return the velocity of the particle relative to c (unitless).
  */
double Material::_beta(double energy_, double mass_){ 
	return std::sqrt(1-pow(mass_/(energy_+mass_), 2.0)); 
}

/** Return the maximum kinetic energy which may be transferred to a free electron per collision.
  * \param[in] beta_ The velocity of the incoming particle relative to c (unitless).
  * \param[in] gamma_ The Lorentz factor for the incoming particle (unitless).
  * \param[in] mass_ The mass of the incoming particle (MeV/c^2).
  * \return the maximum kinetic energy which may be transferred to a free electron (MeV).
  */
double Material::_tmax(double beta_, double gamma_, double mass_){ 
	return (2.0*electron_RME*std::pow(beta_*gamma_, 2.0) / (1.0 + (2.0*gamma_*electron_RME/mass_) + std::pow(electron_RME/mass_, 2.0))); 
}

/** Return the energy for a particle.
  * \param[in] beta2_ The squared velocity of the incoming particle relative to c (unitless).
  * \param[in] mass_ The mass of the incoming particle (MeV/c^2).
  * \return the energy of the particle (MeV).
  */
double Material::_energy(double beta2_, double mass_){ 
	return mass_*((1.0/std::sqrt(1.0-beta2_))-1.0); 
}

/** Return the mean ionization potential.
  * See S. Seltzer, International Journal of Applied Radiation and Isotopes 33, 1189 (1982).
  * \param[in] Z_ The atomic number of the element of interest.
  * \return the ionization potential (MeV) or -1 if the atomic charge is not in the range [1,100].
  */
double Material::_ionpot(unsigned int Z_){ 
	if(Z_ == 1){ return 19.2E-6; }
	else if(Z_ == 6){ return 81E-6; }
	else if(Z_ == 7){ return 82E-6; }
	else if(Z_ == 8){ return 106E-6; }
	else if(Z_ == 9){ return 112E-6; }
	else if(Z_ == 17){ return 180E-6; }
	return (Z_ > 0 && Z_ <= 100 ? potentials[Z_-1]*1.13E-6 : -1.0);
}

/** Return the density effect correction term.
  * See R. Sternheimer, Phys. Rev. B 3, 3681 (1971).
  * \param[in] energy_ The energy of the incident particle in MeV.
  * \return the density effect correction (MeV).
  */
double Material::_density(double energy_){
	return 0.0;
}

/** Default constructor.
  */
Material::Material(){
	_initialize();
}

/** Constructor which initializes all data arrays to the correct size.
  * \param[in] num_elements_ The number of unique elements in this material.
  */
Material::Material(unsigned int num_elements_){
	_initialize();
	Init(num_elements_);
}

/** Destructor.
  */
Material::~Material(){
	if(init){
		delete[] num_per_molecule;
		delete[] element_Z;
		delete[] element_A;
		delete[] element_I;
		delete[] weight;
	}
}

/** Initialize all element variable arrays.
  * \param[in] num_elements_ The number of unique elements in this material.
  * \return true upon success or false if the object has already been initialized.
  */
bool Material::Init(unsigned int num_elements_){
	if(init){ return false; }
	num_elements = num_elements_;
	num_per_molecule = new unsigned int[num_elements];
	element_Z = new double[num_elements];
	element_A = new double[num_elements];
	element_I = new double[num_elements];
	weight = new double[num_elements];
	for(unsigned int i = 0; i < num_elements; i++){
		num_per_molecule[i] = 0;
		element_Z[i] = 0.0;
		element_A[i] = 0.0;
		element_I[i] = 0.0;
		weight[i] = 0.0;
	}

	rad_length = 0.0;
		
	init = true;
	return true;
}

/** Set the Z and A of all of the unique elements in this material. All input arrays should be the same size.
  * \param[in] num_per_molecule_ Pointer to an array containing the number of each element in one molecule of this material.
  * \param[in] Z_ Pointer to an array containing the atomic charge of each element in one molecule of this material.
  * \param[in] A_ Pointer to an array containing the molar mass (in g/mol) of each element in one molecule of this material.
  */
void Material::SetElements(unsigned int *num_per_molecule_, double *Z_, double *A_){
	// Set all unique elements.
	for(unsigned int i = 0; i < num_elements; i++){
		num_per_molecule[i] = num_per_molecule_[i];
		total_elements += num_per_molecule_[i];
		element_Z[i] = Z_[i];
		element_A[i] = A_[i];
		element_I[i] = _ionpot(element_Z[i]);
	}

	_calculate();
}

/** Load a material from a file.
  * \param[in] filename_ The name of the material file to load.
  * \return true upon success or false if loading fails for whatever reason.
  */
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
				element_I[i] = _ionpot(element_Z[i]);
				
				total_elements += num_per_molecule[i];
			}
			
			_calculate();
		}
		count++;
	}
	if(count <= 2){ return false; }

	return true;
}

/** Compute and return the stopping power for an incoming particle in this material.
  * See C. Amsler et al., PL B667, 1 (2008). 
  * \param[in] energy_ The energy of the incoming particle (in MeV).
  * \param[in] Z_ The atomic charge of the incoming particle.
  * \param[in] mass_ The mass of the incoming particle (in MeV/c^2).
  * \return the stopping power (MeV/m).
  */
double Material::StopPower(double energy_, double Z_, double mass_){
	double beta = _beta(energy_, mass_);
	double gamma = 1.0 / std::sqrt(1.0 - beta*beta);
	double tmax = _tmax(beta, gamma, mass_);
	
	double output = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){
		output += weight[i] * (coeff * std::pow(Z_/beta, 2.0) * element_Z[i] / element_A[i]) * (0.5 * std::log(2.0*electron_RME*std::pow(beta*gamma, 2.0)*tmax) - lnIbar - beta*beta - _density(energy_)/2.0);
	}
	
	return output;
}

/** Compute and return the range of an incoming particle in this material.
  * \param[in] energy_ The energy of the incoming particle (in MeV).
  * \param[in] Z_ The atomic charge of the incoming particle.
  * \param[in] mass_ The mass of the incoming particle (in MeV/c^2).
  * \return the range (m).
  */
double Material::Range(double energy_, double Z_, double mass_){
	unsigned int iterations = 100;
	double sum = 0.0;
	double e1, e2;
	double step = energy_/iterations;
	for(unsigned int i = 1; i < iterations; i++){
		e1 = step*i; e2 = step*(i+1);
		sum += 0.5*(1.0/StopPower(e2, Z_, mass_)+1.0/StopPower(e1, Z_, mass_))*step;
	}
	return sum;
}

/** Use Birks' equation to calculate the light output for a particle in this material.
  * \param[in] energy_ The energy of the incoming particle (in MeV).
  * \param[in] Z_ The atomic charge of the incoming particle.
  * \param[in] mass_ The mass of the incoming particle (in MeV/c^2).
  * \param[in] L0_ Light output efficiency of this material (in ??/MeV).
  * \param[in] kB_ Adjustable parameter used for fitting to data (in m/MeV).
  * \param[in] C_ Adjustable parameter used for fitting to data (in m^2/MeV^2).
  * \return the light output (variable unit, takes the unit from the numerator of L0_).
  */
double Material::Birks(double startE_, double energy_, double Z_, double mass_, double L0_, double kB_, double C_/*=0.0*/){
	unsigned int iterations = 100;
	double sum = 0.0;
	double dedx1, dedx2;
	double igrand1, igrand2;
	double step = (energy_-startE_)/iterations;
	dedx1 = -1*StopPower(startE_, Z_, mass_);
	for(unsigned int i = 0; i <= iterations; i++){
		dedx2 = -1*StopPower(startE_+step*(i+1), Z_, mass_);
		igrand1 = 1.0/(1.0 + kB_*dedx1 + C_*dedx1*dedx1);
		igrand2 = 1.0/(1.0 + kB_*dedx2 + C_*dedx2*dedx2);
		sum += 0.5*(igrand1 + igrand2)*step;
		dedx1 = dedx2;
	}
	return L0_*sum;
}

/** Print useful parameters about this material for debugging purposes.
  */ 
void Material::Print(){
	std::cout << " Name: " << vikar_name << std::endl;
	std::cout << "  Number of unique elements: " << num_elements << std::endl;
	for(unsigned int i = 0; i < num_elements; i++){
		std::cout << "   Element " << i+1 << ": N = " << num_per_molecule[i] << ", Z = " << element_Z[i];
		std::cout << ", A = " << element_A[i] << ", I = " << element_I[i] << ", Wn = " << weight[i] << std::endl;
	}
	std::cout << "  Average Z: " << avgZ << "\n";
	std::cout << "  Average A: " << avgA << "\n";
	std::cout << "  Density: " << density << " g/cm^3\n";
	std::cout << "  Molar Mass: " << Mmass << " g/mol\n";
	std::cout << "  Rad Length: " << rad_length << " mg/cm^2\n";
	std::cout << "  Ibar: " << std::exp(lnIbar) << "\n";
}

/** Write useful parameters about this material to a file for debugging purposes.
  * \param[in] file_ Pointer to an output stream file which is open for writing.
  */ 
void Material::Print(std::ofstream *file_){
	if(!file_ || !file_->good()){ return; }
	(*file_) << "Name\t" << vikar_name << std::endl;
	(*file_) << "NumElements\t" << num_elements << std::endl;
	for(unsigned int i = 0; i < num_elements; i++){
		(*file_) << "Element" << i+1 << "\tN=" << num_per_molecule[i] << ", Z=" << element_Z[i];
		(*file_) << ", A=" << element_A[i] << ", I=" << element_I[i] << ", Wn = " << weight[i] << std::endl;
	}
	(*file_) << "AvgZ\t" << avgZ << "\n";
	(*file_) << "AvgA\t" << avgA << "\n";
	(*file_) << "Density\t" << density << " g/cm^3\n";
	(*file_) << "  Molar Mass: " << Mmass << " g/mol\n";
	(*file_) << "  Rad Length: " << rad_length << " mg/cm^2\n";
	(*file_) << "  Ibar: " << std::exp(lnIbar) << "\n";
}

/////////////////////////////////////////////////////////////////////
// Particle
/////////////////////////////////////////////////////////////////////

bool Particle::SetMaterial(Material *mat_, const double &Ebeam_, const double &Espread_){
	if(Z <= 0.0 || init){ return false; }
	mat = mat_;
	maxE = Ebeam_ + 2*Espread_;
	
	if(table.Init(100, 0.1, maxE, Z, mass, mat)){
		init = true;
		return true;
	}
	return false;
}

/// Get the energy of a particle stopped in distance range_.
double Particle::GetTableEnergy(const double &range_){ 
	if(!init){ return -1.0; }
	return table.GetEnergy(range_);
}

/// Get the range of the particle given an initial energy_.
double Particle::GetTableRange(const double &energy_){ 
	if(!init){ return -1.0; }
	return table.GetRange(energy_);
}

/** Get the final energy of a particle with energy_ moving a
  * distance dist_ through the material.
  */
double Particle::GetTableNewE(const double &energy_, const double &dist_, double &dist_traveled){
	if(!init){ return -1.0; }
	return table.GetNewE(energy_, dist_, dist_traveled);
}

/////////////////////////////////////////////////////////////////////
// Target
/////////////////////////////////////////////////////////////////////

Target::Target() : Particle() {
	thickness = 0.0;
	Zthickness = 0.0;
	density = 1.0;
	rad_length = 0.0;
	angle = 0.0;
	physical = new Primitive();
}

Target::Target(unsigned int num_elements_) : Particle() {
	thickness = 0.0;
	Zthickness = 0.0;
	density = 1.0;
	rad_length = 0.0;
	angle = 0.0;
	physical = new Primitive();
}

void Target::SetThickness(double thickness_){ 
	thickness = thickness_;	
	physical->SetSize(1.0, 1.0, thickness/(density*1E5)); // Only the thickness matters
}

void Target::SetRealThickness(double thickness_){ 
	thickness = thickness_*density*1E3;	
	physical->SetSize(1.0, 1.0, thickness_/100.0); // Only the thickness matters
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

/** Get the depth into the target at which the reaction occurs
  * \\param[in] offset_ is the global position where the beam particle originates
  * \\param[in] direction_ is the direction of the beam particle entering the target
  * \\param[out] intersect is the global position where the beam particle intersects the front face of the target
  * \\param[out] interact is the global position where the beam particle reacts inside the target
  */
double Target::GetInteractionDepth(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, Vector3 &interact){
	double t1, t2;
	double zdist = physical->GetApparentThickness(offset_, direction_, intersect, t1, t2); // The target thickness the ray sees
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
