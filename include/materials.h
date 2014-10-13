// materials.h
// Cory Thornsberry

#ifndef MATERIALS_H
#define MATERIALS_H

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

extern const double e_charge, e_mass, e_radius;
extern const double bethe_coefficient;

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
// Target
/////////////////////////////////////////////////////////////////////

class Target{
  private:
  	unsigned int num_elements; // Number of unique elements per molecule of the target
  	unsigned int total_elements; // Total number of elements per molecule in the target
  	unsigned int *num_per_molecule; // Number of each element per molecule
  	double *mean_excitations; // Mean excitations (MeV) for each element
  	double *element_Z, *element_A; // Atomic numbers and atomic masses (u) for each element
  	double Z, A; // Z and A of the target isotope of interest
  	double avgZ, avgA; // Average Z and A of the elements within the target molecule
	double density; // Density of the target (g/cm^3)
	double thickness; // Thickness of the target (mg/cm^2)
	double Zthickness; // Thickness of the target in the z-direction (mg/cm^2)
	double angle; // Angle of target wrt the beam axis (rad)
	double rad_length; // The radiation length of the target (mg/cm^2)
	bool init;
	
	void _initialize();
	
  public:
	Target();
	Target(unsigned int);
	~Target();
	
	bool Init(unsigned int);
	
	void SetZ(double Z_){ Z = Z_; }
	void SetA(double A_){ A = A_; }
	void SetDensity(double density_){ density = density_; }
	void SetThickness(double thickness_){ thickness = thickness_;	}
	void SetAngle(double angle_){
		angle = angle_;	
		Zthickness = dabs(thickness / std::cos(angle));
	}
	void SetElements(unsigned int*, double*, double*);
	
	bool IsInit(){ return init; }
	
	double GetZ(){ return Z; } // Return the Z of the target element of interest
	double GetA(){ return A; } // Return the A of the target element of interest
	double GetAverageZ(){ return avgZ; } // Return the average atomic number of the elements in the molecule
	double GetAverageA(){ return avgA; } // Return the average atomic mass of the elements in the molecule
	double GetRadLength(){ return rad_length; } // Return the radiation length of the target
	double GetDensity(){ return density; } // Return the density of the target (g/cm^3)
	double GetThickness(){ return thickness; } // Return the thickness of the target (mg/cm^2)
	double GetZthickness(){ return Zthickness; } // Return the thickness the beam sees (mg/cm^2)
	double GetRealThickness(){ return thickness/(density*1000.0); } // Return the physical thickness of the target (cm)
	double GetRealZthickness(){ return Zthickness/(density*1000.0); } // Return the physical thickness the beam sees (cm)
	double GetAngle(){ return angle; } // Return the angle of the target wrt the beam axis

	unsigned int GetTotalElements(){ return total_elements; } // Return the total number of elements in the target molecule
	unsigned int GetNumElements(){ return num_elements; } // Return the number of unique elements per target molecule
	
	// Get the real world coordinates of the interaction point within the target
	// Return the Z distance into the target
	double GetInteractionPoint(double, double, Vector3&);
	
	// Return the stopping power of the target
	double GetStoppingPower(double, double, double);
	
	// Return the adjusted energy for a charged particle traversing the target
	double GetEnergy(double, double, double, double);
	
	// Return the maximum range for a charged particle traversing the target
	double GetRange(double, double, double);
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
	bool Init(unsigned int, double, double, double, double, double, double, double, Target*); // Initialize arrays and fill them using ncdedx
	bool UseTable(){ return use_table; }
	bool Set(unsigned int, double, double); // Manually set a data point with an energy and a range
	unsigned int GetEntries(){ return num_entries; } // Return the number of entries in the array
	double GetRange(double); // Get the particle range at a given energy using linear interpolation
	double GetEnergy(double); // Get the particle energy at a given range usign linear interpolation
};

#endif
