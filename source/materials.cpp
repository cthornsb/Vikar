// materials.cpp
// Cory Thornsberry
// Wed Oct 8, 2014

#include "vikar_core.h"
#include "materials.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals 
/////////////////////////////////////////////////////////////////////

const double e_charge = 1.602176565E-19; // C
const double e_mass = 0.510998928; // MeV/c^2
const double e_radius = 2.8179E-15; // m (classical)
const double bethe_coefficient = 4*pi*avagadro*pow(e_radius, 2.0)*e_mass; // MeV*m^2/g

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

bool RangeTable::Init(unsigned int num_entries_, double startE, double stopE, double tgt_dens, double targA, double targZ, double A, double Z){
	if(!_initialize(num_entries_)){ return false; }
	
	double dummy1, dummy2, rangemg;
	
	// Use ncdedx to fill the arrays
	double step = (stopE-startE)/(num_entries-1);
	for(unsigned int i = 0; i < num_entries; i++){
		ncdedx(tgt_dens, targA, targZ, A, Z, (startE+i*step), dummy1, dummy2, rangemg);
		energy[i] = startE + i*step; // Energy in MeV
		range[i] = rangemg/(tgt_dens*1E5); // Range in m
		std::cout << energy[i] << "\t" << range[i] << std::endl;
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
	num_elements = 0;
	total_elements = 0;
	num_per_molecule = NULL;
	element_Z = NULL;
	element_A = NULL;
	init = false;
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
	}
}

bool Material::Init(unsigned int num_elements_){
	if(init){ return false; }
	num_elements = num_elements_;
	num_per_molecule = new unsigned int[num_elements];
	element_Z = new double[num_elements];
	element_A = new double[num_elements];
	for(unsigned int i = 0; i < num_elements; i++){
		num_per_molecule[i] = 0;
		element_Z[i] = 0.0;
		element_A[i] = 0.0;
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
	}
	
	avgZ = 0.0;
	avgA = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){ 
		avgZ += element_Z[i]*num_per_molecule[i];
		avgA += element_A[i]*num_per_molecule[i]; 
	}
	avgZ = avgZ / total_elements;
	avgA = avgA / total_elements;

	// Calculate the radiation length for the target material
	// for use in angular straggling calculations later.
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
	rad_length = 0.0;
	for(unsigned int i = 0; i < num_elements; i++){ 
		rad_length += ((element_A[i]*num_per_molecule[i]/avgA) / radlength(element_A[i],element_Z[i])); 
	}

	rad_length = 1.0/rad_length; 
}

bool Material::ReadMatFile(const char* filename_){
	if(init){ return false; }
	std::ifstream input_file(filename_);
	if(!input_file.good()){ return false; }

	unsigned int num_elements;
	unsigned int *num_per_molecule;
	double *element_Z;
	double *element_A;
	
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
			num_per_molecule = new unsigned int[num_elements];
			element_Z = new double[num_elements];
			element_A = new double[num_elements];
			
			for(unsigned int i = 0; i < num_elements; i++){
				getline(input_file, line); line = Parse(line); element_Z[i] = (double)atof(line.c_str()); 
				getline(input_file, line); line = Parse(line); element_A[i] = (double)atof(line.c_str());  
				getline(input_file, line); line = Parse(line); num_per_molecule[i] = (unsigned int)atoi(line.c_str()); 
			}
		}
		count++;
	}
	if(count <= 2){ return false; }

	Init(num_elements);
	SetElements(num_per_molecule, element_Z, element_A);

	return true;
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
