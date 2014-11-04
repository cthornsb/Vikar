// materials.h
// Cory Thornsberry

#ifndef MATERIALS_H
#define MATERIALS_H

#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

// Constants for use in stopping power calcuations
extern const double electron_RME, proton_RME, neutron_RME;
extern const double bohr_e_radius, e, h, coeff;
extern const double alpha, avagadro;

// Constant coefficients for the shell correction term
extern const double a1, a2, b1, b2, c1, c2;

// Ionization potentials
extern const double ionpot[13];

/////////////////////////////////////////////////////////////////////
// Class declarations
/////////////////////////////////////////////////////////////////////

class Planar;

/////////////////////////////////////////////////////////////////////
// Efficiency
/////////////////////////////////////////////////////////////////////

class Efficiency{
  private:
	double *small_energy, *small_efficiency;
	double *med_energy, *med_efficiency;
	double *large_energy, *large_efficiency;
	unsigned int NsmallEff, NmedEff, NlargeEff;
	bool init_small, init_med, init_large;
	
	bool _read_eff_file(const char* fname, std::vector<double> &energy, std::vector<double> &efficiency);
	
  public:
	Efficiency();
	~Efficiency();
	
	unsigned int GetNsmall(){ return NsmallEff; }
	unsigned int GetNmedium(){ return NmedEff; }
	unsigned int GetNlarge(){ return NlargeEff; }
	
	bool IsSmallInit(){ return init_small; }
	bool IsMediumInit(){ return init_med; }
	bool IsLargeInit(){ return init_large; }
	
	unsigned int ReadSmall(const char*);
	unsigned int ReadMedium(const char*);
	unsigned int ReadLarge(const char*);
	double GetSmallEfficiency(double);
	double GetMediumEfficiency(double);
	double GetLargeEfficiency(double);
};

/////////////////////////////////////////////////////////////////////
// Material
/////////////////////////////////////////////////////////////////////

class Material{
  protected:
	std::string vikar_name; // The name vikar uses to search for this material
  	unsigned int num_elements; // Number of unique elements per molecule of the material
  	unsigned int total_elements; // Total number of elements per molecule in the material
  	unsigned int *num_per_molecule; // Number of each element per molecule
  	double *element_Z, *element_A; // Atomic numbers and atomic masses (u) for each element
  	double *element_I; // Ionization potentials (MeV) 
  	double avgZ, avgA; // Average Z and A of the elements within the material molecule
	double density; // Density of the material (g/cm^3)
	double edens; // Electron density (1/m^3)
	double rad_length; // The radiation length of the material (mg/cm^2)
	double lnIbar; // The natural log of the average ionization potential
	bool init;
	
	void _initialize();

	void _calculate(); // Calculate the average atomic charge, mass, and ionization potential for the material

  public:
  	Material();
  	Material(unsigned int);
	~Material();
	
	bool Init(unsigned int);
	
	void SetDensity(double density_){ density = density_; }
	void SetElements(unsigned int*, double*, double*);
	void SetName(std::string name_){ vikar_name = name_; }

	bool IsInit(){ return init; }
	
	double GetAverageZ(){ return avgZ; } // Return the average atomic number of the elements in the molecule
	double GetAverageA(){ return avgA; } // Return the average atomic mass of the elements in the molecule
	double GetRadLength(){ return rad_length; } // Return the radiation length of the material (mg/cm^2)
	double GetDensity(){ return density; } // Return the density of the material (g/cm^3)
	std::string GetName(){ return vikar_name; }

	unsigned int GetTotalElements(){ return total_elements; } // Return the total number of elements in the material molecule
	unsigned int GetNumElements(){ return num_elements; } // Return the number of unique elements per material molecule
	
	// Read a material file
	bool ReadMatFile(const char* filename_);
	
	double GetEdensity(){ return edens; } // Return the electron density (1/m^3) for the material
	double GetLNibar(){ return lnIbar; } // Return the natural log of the average ionization potential
	
	// Return the effective atomic charge (unitless) for a given value of beta_ and Z_
	double Zeff(double beta_, double Z_){ return Z_*(1-std::exp(-125.0*beta_/pow(Z_, 2.0/3.0))); }

	// Return the value of beta^2 (unitless) for a particle with a given energy_ and mass_
	// energy_ in MeV and mass_ in MeV/c^2
	double Beta2(double energy_, double mass_){ return std::sqrt(2*energy_/mass_); }

	// Return the ionization potential (MeV) for a particle with a given atomic charge Z_
	// Z_ is the atomic number of the ion of interest
	double GetIonPot(unsigned int Z_);

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
	
	void Print();
};

/////////////////////////////////////////////////////////////////////
// Target
/////////////////////////////////////////////////////////////////////

class Target{
  private:
	double thickness; // Thickness of the target (mg/cm^2)
	double Zthickness; // Thickness of the target in the z-direction (mg/cm^2)
	double density; // Density of the target (g/cm^3)
	double rad_length; // The radiation length of the material (mg/cm^2)
	double angle; // Angle of target wrt the beam axis (rad)
  	double Z, A; // Z and A of the material isotope of interest
	
	Planar *physical; // The physical target geometry
	
  public:
	Target();
	Target(unsigned int);

	void SetZ(double Z_){ Z = Z_; }
	void SetA(double A_){ A = A_; }		
	void SetThickness(double thickness_);
	void SetAngle(double angle_);
	void SetDensity(double density_){ density = density_; }
	void SetRadLength(double rad_length_){ rad_length = rad_length_; }
	
	double GetZ(){ return Z; } // Return the Z of the material element of interest
	double GetA(){ return A; } // Return the A of the material element of interest
	double GetThickness(){ return thickness; } // Return the thickness of the target (mg/cm^2)
	double GetZthickness(){ return Zthickness; } // Return the thickness the beam sees (mg/cm^2)
	double GetRealThickness(){ return thickness/(density*1E5); } // Return the physical thickness of the target (m)
	double GetRealZthickness(){ return Zthickness/(density*1E5); } // Return the physical thickness the beam sees (m)
	double GetAngle(){ return angle; } // Return the angle of the target wrt the beam axis
	double GetDensity(){ return density; }
	double GetRadLength(){ return rad_length; }
	Planar *GetPlanar(){ return physical; } // Return a pointer to the 3d geometry object
	
	// Get the depth into the target at which the reaction occurs
	// offset_ is the global position where the beam particle originates
	// direction_ is the direction of the beam particle entering the target
	// intersect is the global position where the beam particle intersects the front face of the target
	// interact is the global position where the beam particle reacts inside the target
	double GetInteractionDepth(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, Vector3 &interact);

	// Determine the new direction of a particle inside the target due to angular straggling
	bool AngleStraggling(const Vector3 &direction_, double A_, double Z, double E_, Vector3 &new_direction);
};

/////////////////////////////////////////////////////////////////////
// RangeTable
/////////////////////////////////////////////////////////////////////

class RangeTable{
  private:
	double *energy, *range;
	unsigned int num_entries;
	bool use_table;
	
	bool _initialize(unsigned int);
	
  public:
	RangeTable(){ use_table = false; }
	RangeTable(unsigned int);
	~RangeTable();
	
	bool Init(unsigned int); // Initialize arrays but do not fill them
	bool Init(unsigned int num_entries_, double startE, double stopE, double A, double Z, Material *mat); // Initialize arrays and fill them using Material
	bool UseTable(){ return use_table; }
	bool Set(unsigned int, double, double); // Manually set a data point with an energy and a range
	unsigned int GetEntries(){ return num_entries; } // Return the number of entries in the array
	double GetRange(double); // Get the particle range at a given energy using linear interpolation
	double GetEnergy(double); // Get the particle energy at a given range usign linear interpolation
};

#endif
