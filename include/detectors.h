// vikar_lib.h
// Cory Thornsberry

#ifndef DETECTORS_H
#define DETECTORS_H

#include <string>
#include <cmath>
#include <fstream>
#include <vector>

#include "vikar_core.h"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

extern const Vector3 zero_vector;

/////////////////////////////////////////////////////////////////////
// Class declarations
/////////////////////////////////////////////////////////////////////

class Material;

/////////////////////////////////////////////////////////////////////
// NewVIKARDet
/////////////////////////////////////////////////////////////////////

struct NewVIKARDet{
	float data[9];
	std::string type;
	std::string subtype;
	std::string material;
	
	NewVIKARDet(){ 
		for(unsigned int i = 0; i < 9; i++){ data[i] = 0.0; }
		type = "unknown";
		subtype = "unknown"; 
		material = "none";
	}
	NewVIKARDet(std::string input_){
		for(unsigned int i = 0; i < 9; i++){ data[i] = 0.0; }
		SetValues(input_);
	}
	
	void SetValues(std::string input){
		std::string temp_str = "";
		unsigned int current_index = 0;
		bool done = false;
		for(unsigned int i = 0; i < input.size(); i++){
			if(input[i] == '#'){ break; } // No data should come after a comment marker
			else if(i == input.size()-1){
				temp_str += input[i];
				done = true;
			}
			else if((input[i] == ' ' || input[i] == '\t' || input[i] == '\n') && temp_str != ""){ done = true; }
			else{ temp_str += input[i]; }
			
			if(done){
				if(current_index <= 5){ data[current_index] = atof(temp_str.c_str()); }
				else if(current_index == 6){ type = temp_str; }
				else if(current_index == 7){ subtype = temp_str; }
				else if(current_index == 8){ data[6] = atof(temp_str.c_str()); }
				else if(current_index == 9){ data[7] = atof(temp_str.c_str()); }
				else if(current_index == 10){ data[8] = atof(temp_str.c_str()); }
				else if(current_index == 11){ material = temp_str; }
				else{ break; }
				current_index++;
				temp_str = "";
				done = false;
			}
		}
	}
	
	std::string DumpDet(){
		std::stringstream stream;
		stream << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5];
		stream << "\t" << type << "\t" << subtype << "\t" << data[6] << "\t" << data[7] << "\t" << data[8];
		return stream.str();
	}
};

/////////////////////////////////////////////////////////////////////
// Planar
/////////////////////////////////////////////////////////////////////

class Planar{  
    private:
	bool need_set;

	// Set the global face coordinates (wrt global origin)
	// Each vertex is the center coordinate of one of the faces 	
	void _set_face_coords();
	 	
    protected:
    bool use_material;
    unsigned int material_id;
	Vector3 position; // Center, cartesian position
	Vector3 detX, detY, detZ; // Local face unit vectors
	Vector3 GlobalFace[6]; // 6 global face coordinates
	double length, width, depth; // Physical size (meters)
	double theta, phi, psi; // Local rotation (radians)
	bool small, med, large;
	bool use_recoil;
	bool use_eject;
	std::string type, subtype;
    
    public:
	Planar();
	
	unsigned int GetMaterial(){ return material_id; }

	// Return the unit vector of one of the faces. Returns the zero vector if the face is undefined
	void GetUnitVector(unsigned int face_, Vector3 &unit);

	void GetPosition(Vector3 &pos){ pos = position; }

	double GetX(){ return position.axis[0]; }

	double GetY(){ return position.axis[1]; }

	double GetZ(){ return position.axis[2]; }

	double GetWidth(){ return width; }

	double GetLength(){ return length; }

	double GetDepth(){ return depth; }

	std::string GetType(){ return type; }

	std::string GetSubtype(){ return subtype; }

	// Return the local bar frame coordinates of a global coordinate
	void GetLocalCoords(const Vector3&, double&, double&, double&);
	
	bool IsSmall(){ return small; }

	bool IsMedium(){ return med; }

	bool IsLarge(){ return large; }

	bool IsRecoilDet(){ return use_recoil; }

	bool IsEjectileDet(){ return use_eject; }

	bool UseMaterial(){ return use_material; }
	
	void SetMaterial(unsigned int material_id_){ 
		material_id = material_id_; 
		use_material = true;
	}

	void SetType(std::string type_){ type = type_; }

	void SetSubtype(std::string subtype_){ subtype = subtype; }

	void SetSmall(){ length = 0.60; width = 0.03; depth = 0.03; small = true; med = false; large = false; need_set = true; }

	void SetMedium(){ length = 1.20; width = 0.05; depth = 0.03; small = false; med = true; large = false; need_set = true; }

	void SetLarge(){ length = 2.00; width = 0.05; depth = 0.05; small = false; med = false; large = true; need_set = true; }

	// Set the physical size of the bar
	// For known bar sizes, it is better to use SetSmall/SetMedium/SetLarge methods
	// Unknown bar types will not include efficiency data
	void SetSize(double, double, double);
	
	// Set the position of the center of the bar using a 3d vector (in meters)
	void SetPosition(const Vector3&);
	
	// Set the cartesian position of the center of the bar in 3d space (in meters)
	void SetPosition(double, double, double);

	// Set the polar position of the center of the bar in 3d space (meters and radians)
	void SetPolarPosition(double, double, double);
	
	// Get the local unit vectors using 3d matrix rotation
	// X and Y are face axes, Z is the axis into or out of the bar
	void SetRotation(double, double, double);

	// Manually set the local bar unit vectors
	// I would advise against this, use SetRotation instead
	// Note: This method does not update the bar rotation values
	//  and should therefore only be used for testing and debugging
	void SetUnitVectors(const Vector3&, const Vector3&, const Vector3&);
	
	void SetRecoil(bool input_=true){ use_recoil = input_; }
	
	void SetEjectile(bool input_=true){ use_eject = input_; }
	
	// Check if a point (in local coordinates) is within the bounds of the primitive
	// Return true if the coordinates are within the primitive and false otherwise
	bool CheckBounds(unsigned int face_, double x_, double y_, double z_);

	// Find if a ray (from the origin) intersects the infinite
	// plane defined by one of the faces of the VANDLE bar
	// face 0 is along the +z local axis
	// face 1 is along the +x local axis
	// face 2 is along the -z local axis
	// face 3 is along the -x local axis
	// face 4 is along the +y local axis
	// face 5 is along the -y local axis	
	bool PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, Vector3 &P);
	
	// Determine if a ray from the origin intersected a face of the bar
	// and return the face id number. Return -1 if no intersection was found
	int FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz);
	
	// Check for an intersect with the top or bottom faces.
	// In order to save time, you would normally not want to check
	// this because the PMTs will block the top and bottom faces
	int TopBotIntersect(const Vector3 &offset_, const Vector3 &direction, Vector3 &intersect, double &px, double &py, double &pz);
	
	// Calculate the intersection of a ray of the form (offset_ + t * direction_) with this primitive shape
	// offset_ is the point where the ray originates wrt the global origin
	// direction_ is the direction of the ray wrt the global origin
	// P1 is the first intersection point in global coordinates
	// P2 is the second intersection point in global coordinates
	// Return true if the primitive is intersected, and false otherwise
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &P2, int &face1, int &face2, double &px, double &py, double &pz);
	
	// Trace a ray through the detector and calculate the thickness it sees between two faces (f1 and f2)
	// Return -1 if the ray does not travel through both faces
	double GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned int f1_, unsigned int f2_, Vector3 &intersect1, Vector3 &intersect2);
	
	// Dump raw cartesian face vertex data.
	// This returns a string containing the vertex coordinates of the
	// centers of all six faces of the VANDLE bar and its center coordinate
	std::string DumpVertex();
	
	// Dump VIKAR detector format string
	// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
	std::string DumpDet();
};

/////////////////////////////////////////////////////////////////////
// Wall
/////////////////////////////////////////////////////////////////////

class Wall: public Planar{
    private:
	unsigned int num_bars;
	Planar *bars;
	bool init; 
	
    public:
	Wall(){ bars = NULL; init = false; num_bars = 0; }
	~Wall(){ if(init){ delete[] bars; } }
	
	// Initialize the wall object.
	// Spacing is the physical distance between the edges of bars.
	// This method leaves a distance of spacing/2 on each side of the wall
	// so that walls having the same bar spacing will align properly
	void Initialize(unsigned int, double, double, double, double);
	
	// Return the total number of detectors
	unsigned int GetNumBars(){ return num_bars; }
	
	// Return a pointer to a bar in this wall
	// Return null if the bar does not exist
	Planar *GetBar(unsigned int);
	
	// Return a pointer to the array of planar detectors
	Planar *GetBars(){ return bars; }
	
	// Test the planar wall class
	void WallTest();
	
	// Dump VIKAR detector format string
	// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
	std::string DumpDetWall();
};

/// Read in a detector efficiency file
unsigned int ReadEffFile(const char*, double*, double*);

// Read NewVIKAR detector file and load bars into an array
// Returns the number of detectors loaded from the file
// Assumes the following detector file format for each bar in file
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
unsigned int ReadDetFile(const char* fname_, std::vector<Planar*> &bar_vector);

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates one output root file named 'mcarlo.root'
// fwhm_ (m) allows the use of a gaussian particle "source". If fwhm_ == 0.0, a point source is used
// angle_ (rad) allows the rotation of the particle source about the y-axis
unsigned int TestDetSetup(Planar *bar_array, unsigned int num_bars, unsigned int num_trials, bool WriteRXN_, double fwhm_=0.0, double angle_=0.0);

#endif
