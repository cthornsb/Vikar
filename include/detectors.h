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
	bool is_cylinder; // TEMPORARY!!!
	bool use_recoil;
	std::string type, subtype;
    
    public:
	Planar();
	
	unsigned int GetMaterial(){ return material_id; }
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
	void GetLocalCoords(const Vector3&, double&, double&, double&);
	
	bool IsCylinder(){ return is_cylinder; }
	bool IsSmall(){ return small; }
	bool IsMedium(){ return med; }
	bool IsLarge(){ return large; }
	bool IsRecoilDet(){ return use_recoil; }
	bool UseMaterial(){ return use_material; }
	
	void SetMaterial(unsigned int material_id_){ material_id = material_id_; }
	void SetCylinder(){ is_cylinder = true; }
	void SetType(std::string type_){ type = type_; }
	void SetSubtype(std::string subtype_){ subtype = subtype; }
	void SetSmall(){ length = 0.60; width = 0.03; depth = 0.03; small = true; med = false; large = false; need_set = true; }
	void SetMedium(){ length = 1.20; width = 0.05; depth = 0.03; small = false; med = true; large = false; need_set = true; }
	void SetLarge(){ length = 2.00; width = 0.05; depth = 0.05; small = false; med = false; large = true; need_set = true; }
	void SetSize(double, double, double);
	void SetPosition(const Vector3&);
	void SetPosition(double, double, double);
	void SetPolarPosition(double, double, double);
	void SetRotation(double, double, double);
	void SetUnitVectors(const Vector3&, const Vector3&, const Vector3&);
	void SetRecoil(bool input_=true){ use_recoil = input_; }
	
	bool CheckBounds(unsigned int face_, double x_, double y_, double z_);
	bool PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, Vector3 &P);
	int FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz);
	int TopBotIntersect(const Vector3 &offset_, const Vector3 &direction, Vector3 &intersect, double &px, double &py, double &pz);
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &P2, int &face1, int &face2, double &px, double &py, double &pz);
	double GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned int f1_, unsigned int f2_, Vector3 &intersect1, Vector3 &intersect2);
	std::string DumpVertex();
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
	void Initialize(unsigned int, double, double, double, double);
	unsigned int GetNumBars(){ return num_bars; }
	Planar *GetBar(unsigned int);
	Planar *GetBars(){ return bars; }
	void WallTest();
	std::string DumpDetWall();
};

unsigned int ReadEffFile(const char*, double*, double*);
unsigned int ReadDetFile(const char* fname_, std::vector<Planar*> &bar_vector);
unsigned int TestDetSetup(Planar *bar_array, unsigned int num_bars, unsigned int num_trials);

#endif
