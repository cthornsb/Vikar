// materials.h
// Cory Thornsberry

#ifndef MATERIALS_H
#define MATERIALS_H

#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

// Constants for use in stopping power calcuations
extern const double electron_RME;
extern const double proton_RME;
extern const double neutron_RME;
extern const double bohr_e_radius;
extern const double e_charge;
extern const double alpha, avagadro;
extern const double mev2amu;
extern const double mev2kg;

// Ionization potentials
extern const double ionpot[13];

/////////////////////////////////////////////////////////////////////
// Class declarations
/////////////////////////////////////////////////////////////////////

class Primitive;
class Efficiency;
class Material;
class Particle;
class Target;
class RangeTable;

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
	std::string vikar_name; /// The name vikar uses to search for this material.
  	unsigned int num_elements; /// Number of unique elements per molecule of the material.
  	unsigned int total_elements; /// Total number of elements per molecule in the material.
  	unsigned int *num_per_molecule; /// Number of each element per molecule.
  	double *element_Z, *element_A; /// Atomic numbers and atomic masses (u) for each element.
  	double *element_I; /// Ionization potentials (MeV).
  	double *weight; /// The fractional weight of each element.
  	double avgZ, avgA; /// Average Z and A of the elements within the material molecule.
	double density; /// Density of the material (g/cm^3).
	double Mmass; /// Molar mass of material (g/mol).
	double rad_length; /// The radiation length of the material (mg/cm^2).
	double lnIbar; // The natural log of the average ionization potential.
	double coeff; /// The leading coefficient of the Bethe-Bloch equation (MeV * g / (mol * m)).
	bool init; /// Set to true if this material has been initialized correctly.
	bool use_eloss; /// Set to true if this material is able to do dE/dx calculations.

	///Initialize all variables with default values.
	void _initialize();

	///Calculate the average atomic charge, mass, and ionization potential for the material.
	void _calculate();

	///Return the velocity of a particle with a given energy and mass relative to c.
	double _beta(double energy_, double mass_);

	///Return the maximum kinetic energy which may be transferred to a free electron per collision.
	double _tmax(double beta_, double gamma_, double mass_);
	
	///Return the energy for a particle.
	double _energy(double beta2_, double mass_);
	
	///Return the mean ionization potential.
	double _ionpot(unsigned int Z_);

	///Return the density effect correction term.
	double _density(double energy_);

  public:
	/// Default constructor.
  	Material();
  	
	/// Constructor which initializes all data arrays to the correct size.
  	Material(unsigned int);
  	
  	/// Destructor.
	~Material();
	
	/// Initialize all element variable arrays.
	bool Init(unsigned int);
	
	/// Set a flag to either use or not use energy loss calculations.
	void SetUseFlag(bool state_=true){ use_eloss = state_; }
	
	/// Set the density of the material (g/cm^3).
	void SetDensity(double density_){ density = density_; }
	
	/// Set the molar mass of the material (g/mol).
	void SetMolarMass(double Mmass_){ Mmass = Mmass_; }
	
	/// Set the number of unique elements in this material.
	void SetElements(unsigned int*, double*, double*);
	
	/// Set the Vikar name of this material.
	void SetName(std::string name_){ vikar_name = name_; }

	/// Return true if this material is initialized and false otherwise.
	bool IsInit(){ return init; } 
	
	/// Return the average atomic number of the elements in the molecule.
	double GetAverageZ(){ return avgZ; } 
	
	/// Return the average atomic mass of the elements in the molecule.
	double GetAverageA(){ return avgA; } 
	
	/// Return the radiation length of the material (mg/cm^2).
	double GetRadLength(){ return rad_length; } 
	
	/// Return the density of the material (g/cm^3).
	double GetDensity(){ return density; } 
	
	/// Return the natural log of the average ionization potential.
	double GetLNibar(){ return lnIbar; } 
	
	/// Return the molar mass of the material (g/mol).
	double GetMolarMass(){ return Mmass; } 
	
	/// Return the Vikar name of this material.
	std::string GetName(){ return vikar_name; } 

	/// Return the total number of elements in the material molecule.
	unsigned int GetTotalElements(){ return total_elements; } 

	/// Return the number of unique elements per material molecule.
	unsigned int GetNumElements(){ return num_elements; } 
	
	/// Load a material from a file.
	bool ReadMatFile(const char* filename_);
	
	/// Compute and return the stopping power for an incoming particle in this material.
	double StopPower(double energy_, double Z_, double mass_);

	/// Compute and return the range of an incoming particle in this material.
	double Range(double energy_, double Z_, double mass_);
	
	/// Use Birks' equation to calculate the light output for a particle in this material.
	double Birks(double startE_, double energy_, double Z_, double mass_, double L0_, double kB_, double C_=0.0);
	
	/// Print useful parameters about this material for debugging purposes.
	void Print();
	
	/// Write useful parameters about this material to a file for debugging purposes.
	void Print(std::ofstream *file_);
};

/////////////////////////////////////////////////////////////////////
// RangeTable
/////////////////////////////////////////////////////////////////////

class RangeTable{
  private:
	double *energy; /// Array for storing energy values.
	double *dedx; /// Array for storing stopping power.
	double *range; /// Array for storing range values.
	double *birks; /// Array for storing light response.
	double step; /// Energy step size (MeV).
	unsigned int num_entries; /// Number of table array entries.
	bool use_table; /// True if the table is to be used for energy loss calculations.
	bool use_birks; /// True if the birks light response table may be used for calculations.
	
	/// Initialize range table arrays.
	bool _initialize(const unsigned int &num_entries_);
	
	/// Interpolate between two points
	double _interpolate(double *x_, double *y_, const double &val_);
	
  public:
  	/// Default constructor.
	RangeTable();
	
	/// Constructor to set the number of table entries.
	RangeTable(const unsigned int &num_entries_);
	
	/// Destructor.
	~RangeTable();
	
	/// Initialize arrays but do not fill them.
	bool Init(const unsigned int &num_entries_); 

	/// Initialize arrays and fill them using Material.
	bool Init(const unsigned int &num_entries_, const double &startE_, const double &stopE_, const double &Z_, const double &mass_, Material *mat_); 

	/// Initialize the birks light response array.
	bool InitBirks(double L0_, double kB_, double C_=0.0);

	/// Return true if the range table is to be used for energy loss and false otherwise.
	bool UseTable(){ return use_table; }

	/// Manually set a data point with an energy and a range.
	bool Set(const unsigned int &pt_, const double &energy_, const double &range_); 

	/// Return the number of entries in the array.
	unsigned int GetEntries(){ return num_entries; } 

	/// Get the particle range at a given energy using linear interpolation.
	double GetRange(const double &energy_); 

	/// Get the particle energy at a given range using linear interpolation.
	double GetEnergy(const double &range_); 

	/// Get the scintillator light response due to a particle traversing a material with given kinetic energy.
	double GetKEfromLR(const double &energy_);
	
	/// Get the kinetic energy of a particle which produces a given light response in a scintillator.
	double GetLRfromKE(const double &response_);

	/// Get the new energy of a particle traversing a distance through a material.
	double GetNewE(const double &energy_, const double &dist_); 

	/// Get the new energy of a particle traversing a distance through a material.
	double GetNewE(const double &energy_, const double &dist_, double &dist_traveled);
	
	/// Return the range and energy for an entry in the table.
	bool GetEntry(const unsigned int &entry_, double &E, double &R);
	
	/// Print range table entries to the screen.
	void Print();
};

/////////////////////////////////////////////////////////////////////
// Particle
/////////////////////////////////////////////////////////////////////

class Particle{
  protected:
	double A, Z; /// Mass and charge number of the particle
	double maxE; /// Maximum energy in the range table
	double mass; /// Rest mass energy of the particle (MeV/c^2)
	std::string name; /// The name of the particle.
	bool init; /// True if the material and range table have been initialized.
	
	RangeTable table; /// The range table to use for energy loss calculations.
	Material *mat; /// The material to use for energy loss calculations.
	
  public:
	Particle(){ 
		A = 0.0; Z = 0.0; maxE = 0.0; mass = 0.0;
		name = "unknown"; init = false; 
	}
	Particle(const std::string &name_, const double &Z_, const double &A_, const double &BE_A_=0.0){ 
		SetParticle(name_, Z_, A_, BE_A_); 
		init = false;
	}
	
	/**Set the mass number of the particle. This does NOT set 
	  * the actual mass of the particle. To do that, call SetMass().
	  */
	void SetA(const double &A_){ A = A_; }
	
	/// Set the charge number of the particle.
	void SetZ(const double &Z_){ Z = Z_; }
	
	/// Auto-set the mass of the particle in units of MeV/c^2. Assumes that Z and A have been set.
	void SetMass(const double &BE_A_=0.0){ mass = Z*proton_RME + (A-Z)*neutron_RME - BE_A_*A; }
	
	/// Manually set the mass of the particle in units of MeV/c^2.
	void SetMassMeV(const double &mass_){ mass = mass_; }
	
	/// Manually set the mass of the particle in units of amu.
	void SetMassAMU(const double &mass_){ mass = mass_/mev2amu; } 
	
	/// Manually set the mass of the particle in units of kg.
	void SetMassKg(const double &mass_){ mass = mass_/mev2kg; }
	
	/// Setup the particle name, charge number, and mass.
	void SetParticle(const std::string &name_, const double &Z_, const double &A_, const double &BE_A_=0.0){
		Z = Z_; A = A_; name = name_; 
		SetMass(BE_A_);
	} 
	
	/// Set the name of the particle.
	void SetName(std::string name_){ name = name_; }
	
	/**Set the material for the particle to use for energy loss calculations.
	  * This also sets up the range table.
	  */
	bool SetMaterial(Material *mat_, const double &Ebeam_, const double &Espread_);
	
	double GetA(){ return A; } /// Return the atomic mass of the particle
	double GetZ(){ return Z; } /// Return the atomic charge of the particle
	double GetN(){ return A-Z; } /// Return the number of neutrons
	double GetMass(){ return mass; } /// Return the rest mass of the particle (MeV/c^2)
	double GetMassKg(){ return mass*mev2kg; } /// Return the rest mass of the particle (kg).
	double GetMassAMU(){ return mass*mev2amu; } /// Return the rest mass of the particle (amu).
	double GetMaxE(){ return maxE; } /// Return the maximum particle energy (MeV)
	double GetCharge(){ return Z*e_charge; } /// Return the electric charge of the particle (C).
	std::string GetName(){ return name; } /// Return the name of the particle.
	Material *GetMaterial(){ return mat; } /// Return a pointer to the energy loss material.
	RangeTable *GetTable(){ return &table; } /// Return a pointer to the range table.
	
	/// Get the relativistic factor of the particle (unitless).
	double GetGamma(const double &velocity_){ return 1.0/std::sqrt(1.0-velocity_*velocity_/(c*c)); }

	/// Get the particle velocity relative to c from the total energy (MeV).
	double GetBetafromTE(const double &energy_){ return GetVfromTE(energy_)/c; }

	/// Get the particle velocity relative to c from the kinetic energy (MeV).
	double GetBetafromKE(const double &energy_){ return GetVfromKE(energy_)/c; }
	
	/// Get the relativistic velocity of the particle (unitless).
	double GetBeta(const double &velocity_){ return velocity_/c; }	
	
	/// Get the kinetic energy from the total energy (MeV).
	double GetKEfromTE(const double &energy_){ return (energy_ - mass); }	
	
	/// Get the total energy from the kinetic energy (MeV).
	double GetTEfromKE(const double &energy_){ return (energy_ + mass); }
	
	/// Get the kinetic energy from the velocity of the particle (MeV).
	double GetKEfromV(const double &velocity_){ return (GetGamma(velocity_) - 1)*mass; }
	
	/// Get the total energy from the velocity of the particle (MeV).
	double GetTEfromV(const double &velocity_){ return (std::sqrt(GetGamma(velocity_)*GetGamma(velocity_)*velocity_*velocity_/(c*c) + 1.0)*mass); }

	/// Get the kinetic energy from the momentum of the particle (MeV/c).
	double GetKEfromP(const double &momentum_){ return (GetTEfromP(momentum_) - mass); }
	
	/// Get the total energy from the velocity of the particle (MeV).
	double GetTEfromP(const double &momentum_){ return std::sqrt(momentum_*momentum_ + mass*mass); }
	
	/// Get the relativistic momentum from the kinetic energy of the particle (MeV/c).
	double GetPfromKE(const double &energy_){ return GetPfromTE(energy_+mass); }	
	
	/// Get the relativistic momentum from the total energy of the particle (MeV/c).
	double GetPfromTE(const double &energy_){ return std::sqrt(energy_*energy_/(c*c) - mass*mass); }

	/// Get the relativistic momentum from the velocity of the particle (MeV/c).
	double GetPfromV(const double &velocity_){ return GetGamma(velocity_)*mass*velocity_/c; }

	/// Get the relativistic velocity from the kinetic energy of the particle (m/s).
	double GetVfromKE(const double &energy_){ return c*std::sqrt(1.0 - std::pow(1.0/(1.0 + energy_/mass), 2.0)); }

	/// Get the relativistic velocity from the total energy of the particle (m/s).
	double GetVfromTE(const double &energy_){ return c*std::sqrt(1.0 - mass*mass/(energy_*energy_)); }
	
	/// Get the relativistic velocity from the momentum of the particle (MeV/c).
	double GetVfromP(const double &momentum_){ return GetVfromKE(GetKEfromP(momentum_)); }
	
	/// Get the energy of a particle stopped in distance range_.
	double GetTableEnergy(const double &range_);
	
	/// Get the range of the particle given an initial energy_.
	double GetTableRange(const double &energy_);
	
	/**Get the final energy of a particle with energy_ moving a
	  * distance dist_ through the material.
	  */
	double GetTableNewE(const double &energy_, const double &dist_, double &dist_traveled);
};

/////////////////////////////////////////////////////////////////////
// Target
/////////////////////////////////////////////////////////////////////

class Target : public Particle {
  private:
	double thickness; /// Thickness of the target (mg/cm^2).
	double Zthickness; /// Thickness of the target in the z-direction (mg/cm^2).
	double density; /// Density of the target (g/cm^3).
	double rad_length; /// The radiation length of the material (mg/cm^2).
	double angle; /// Angle of target wrt the beam axis (rad).
	
	Primitive *physical; /// The physical target geometry.
	
  public:
  	/// Default constructor.
	Target();
	
	/// Constructor to set the number of target elements.
	Target(unsigned int);

	void SetThickness(double thickness_); /// Set the thickness of the target (in mg/cm^2).
	void SetRealThickness(double thickness_); /// Set the actual thickness of the target (in cm).
	void SetAngle(double angle_); /// Set the angle of the target about the z-axis (in rad).
	void SetDensity(double density_); /// Set the density of the target (in g/cm^3).
	void SetRadLength(double rad_length_){ rad_length = rad_length_; } /// Set the radiation length of the target (in mg/cm^2).
	
	double GetThickness(){ return thickness; } /// Return the thickness of the target (mg/cm^2).
	double GetZthickness(){ return Zthickness; } /// Return the thickness the beam sees (mg/cm^2).
	double GetRealThickness(){ return thickness/(density*1E5); } /// Return the physical thickness of the target (m).
	double GetRealZthickness(){ return Zthickness/(density*1E5); } /// Return the physical thickness the beam sees (m).
	double GetAngle(){ return angle; } /// Return the angle of the target about the z-axis (rad).
	double GetDensity(){ return density; } /// Return the density of the target (g/cm^3).
	double GetRadLength(){ return rad_length; } /// Return the radiation length of the target (mg/cm^2).
	Primitive *GetPrimitive(){ return physical; } /// Return a pointer to the 3d geometry object
	
	/// Get the depth into the target at which the reaction occurs.
	double GetInteractionDepth(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, Vector3 &interact);

	/// Determine the new direction of a particle inside the target due to angular straggling.
	bool AngleStraggling(const Vector3 &direction_, double A_, double Z, double E_, Vector3 &new_direction);
};

#endif
