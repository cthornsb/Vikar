#ifndef VANDMC_HPP
#define VANDMC_HPP

#include <iostream>
#include <string>

// SimpleScan
#include "optionHandler.hpp"

// VANDMC
#include "vandmc_core.hpp"
#include "kindeux.hpp"
#include "materials.hpp"
#include "detectors.hpp"
#include "vandmcStructures.hpp"

///////////////////////////////////////////////////////////////////////////////
// class vandmcParameter
///////////////////////////////////////////////////////////////////////////////

class vandmcParameter{
  public:
	vandmcParameter(){ }

	vandmcParameter(const std::string &name_, const std::string &info_="", const std::string &val_="") : name(name_), info(info_), val(val_) { }

	std::string SetName(const std::string &str){ return (name=str); }

	std::string SetInfo(const std::string &str){ return (info=str); }

	std::string SetValue(const std::string &str){ return (val=str); }

	std::string GetName() const { return name; }

	std::string GetInfo() const { return info; }

	std::string GetValue() const { return val; }

	double ConvertToDouble();
	
	unsigned int ConvertToUlong();
	
	int ConvertToLong();
	
	bool ConvertToBool();
	
	bool Compare(const std::string &str){ return (str==name); }
	
	bool Empty(){ return val.empty(); }

  private:
	std::string name;
	std::string info;
	std::string val;
};

///////////////////////////////////////////////////////////////////////////////
// class vandmcParameterReader
///////////////////////////////////////////////////////////////////////////////

class vandmcParameterReader{
  public:
	vandmcParameterReader(){ initialize(); }
	
	bool Read(const char *fname, bool echo=false);

	bool FindInList(const std::string &str);
	
	size_t FindAllOccurances(const std::string &str, std::vector<vandmcParameter*> &vec);
	
	bool FindDouble(const std::string &str, double &val);
	
	bool FindUlong(const std::string &str, unsigned int &val);
	
	bool FindLong(const std::string &str, int &val);
	
	bool FindBool(const std::string &str, bool &val);
	
	bool FindString(const std::string &str, std::string &val);

  private:
	std::vector<vandmcParameter> validParameters;

	std::vector<vandmcParameter> parameters;

	vandmcParameter *findParam(const std::string &str);

	void initialize();
};

///////////////////////////////////////////////////////////////////////////////
// class vandmc
///////////////////////////////////////////////////////////////////////////////

class vandmc{
  public:
	vandmc(){ initialize(); }

	bool Execute(int argc, char *argv[]);

  private:
	Kindeux kind; // Main kinematics object
	Target targ; // The physical target
	Efficiency bar_eff; // VANDLE bar efficiencies
	std::vector<Primitive*> vandle_bars; // Vector of Primitive detectors

	RangeTable beam_targ; // Range table for beam in target
	RangeTable eject_targ; // Pointer to the range table for ejectile in target
	RangeTable recoil_targ; // Pointer to the range table for recoil in target
	std::vector<RangeTable> eject_tables; // Array of range tables for ejectile in various materials
	std::vector<RangeTable> recoil_tables; // Array of range tables for recoil in various materials

	Particle recoil_part; // Recoil particle
	Particle eject_part;// Ejectile particle
	Particle beam_part; // Beam particle
	
	double hit_x, hit_y, hit_z; // Hit coordinates on the surface of a detector

	Vector3 ZeroVector; // The zero vector
	Vector3 Ejectile, Recoil, Gamma;
	Vector3 HitDetect1, HitDetect2;
	Vector3 RecoilSphere;
	Vector3 EjectSphere;
	Vector3 GammaSphere;
	Vector3 lab_beam_focus; // The focal point for the beam. Non-cylindrical beam particles will originate from this point.
	Vector3 lab_beam_start; // The originating point of the beam particle in the lab frame
	Vector3 lab_beam_trajectory; // The original trajectory of the beam particle before it enters the target
	Vector3 lab_beam_interaction; // The position of the reaction inside the target
	Vector3 lab_beam_stragtraject; // The angular straggled trajectory of the beam particle just before the reaction occurs
	Vector3 targ_surface; // The intersection point between the beam particle and the target surface (wrt beam focus)
	Vector3 interaction; // The interaction point inside the target (wrt beam focus)
	Matrix3 rotation_matrix; // The rotation matrix used to transform vectors from the beam particle frame to the lab frame

	double Zdepth; // Interaction depth inside of the target (m)
	double range_beam;
		
	unsigned int num_materials;
	std::vector<Material> materials; // Array of materials
	
	unsigned int targ_mat_id; // The ID number of the target material
	std::string targ_mat_name; // The name of the target material
	
	unsigned int NRecoilStates;
	std::vector<std::string> AngDist_fname; 
	double *ExRecoilStates;
	double *totXsect;
	double gsQvalue;

	double Ebeam, Ebeam0;
	double ErecoilMod;
	double EejectMod;
	double Egamma;
		
	double beamspot; // Beamspot diameter (m) (on the surface of the target)
	double beamEspread; // Beam energy spread (MeV)
	double beamAngdiv; // Beam angular divergence (radians)

	double timeRes; // Pixie-16 time resolution (s)
	double BeamRate; // Beam rate (1/s)

	unsigned int backgroundRate;
	unsigned int backgroundWait;
	double detWindow;
	bool bgPerDetection;
	
	unsigned int NgoodDetections; // Number of wanted particles detected in ejectile detectors.
	unsigned int Ndetected; // Total number of particles detected in ejectile detectors.
	unsigned int Nwanted; // Number of desired ejectile detections.
	unsigned int Nsimulated; // Total number of simulated particles.
	unsigned int NdetHit; // Total number of particles which collided with a detector.
	unsigned int Nreactions; // Total number of particles which react with the target.
	unsigned int NrecoilHits;
	unsigned int NejectileHits;
	unsigned int NgammaHits;
	unsigned int NvetoEvents;
	int Ndet; // Total number of detectors
	unsigned int NdetRecoil; // Total number of recoil detectors
	unsigned int NdetEject; // Total number of ejectile detectors
	unsigned int NdetGamma; // Total number of gamma detectors
	unsigned int NdetVeto; // Total number of particle vetos
	unsigned int BeamType; // The type of beam to simulate (0=gaussian, 1=cylindrical, 2=halo)
	clock_t timer; // Clock object for calculating time taken and remaining

	bool InverseKinematics;
	bool InCoincidence;
	bool WriteReaction;
	bool NeutronSource;
	bool PerfectDet;
	bool SupplyRates;
	bool BeamFocus;
	bool DoRutherford;
	bool echoMode;
	bool printParams;
	unsigned int ADists;
	
	bool have_recoil_det;
	bool have_ejectile_det;
	bool have_gamma_det;
	bool have_veto_det;

	std::string detector_filename;
	std::string input_filename;
	std::string output_filename;

	vandmcParameterReader reader; // Config file
	optionHandler handler;

	void initialize();

	void titleCard();

	bool readConfig(const char *fname);

	bool setup(int argc, char *argv[]);
	
	void print();
};

#endif
