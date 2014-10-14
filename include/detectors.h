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
// NewVIKARDet
/////////////////////////////////////////////////////////////////////

struct NewVIKARDet{
	float data[9];
	std::string bar_type;
	
	NewVIKARDet(){ bar_type = "unknown"; }
	NewVIKARDet(float data_[9], std::string bar_type_){
		SetValues(data_, bar_type_);
	}
	
	void SetValues(float data_[9], std::string bar_type_){
		bar_type = bar_type_;
		for(unsigned short i = 0; i < 9; i++){ 
			if(bar_type == "small" || bar_type == "medium" || bar_type == "large"){
				if(i < 6){ data[i] = data_[i]; }
				else{ data[i] = 0.0; }
			}
			else{ data[i] = data_[i]; }
		}
	}
};

class Planar{  
    private:
	bool need_set;
	
	void _set_face_coords();
	 	
    protected:
	Vector3 position; // Center, cartesian position
	Vector3 barX, barY, barZ; // Local face unit vectors
	Vector3 GlobalFace[6]; // 6 global face coordinates
	//Vector3 GlobalVertex[6][4]; // Vertex data for the faces
	double length, width, depth; // Physical bar size (meters)
	double theta, phi, psi; // Local bar rotation (radians)
	bool small, med, large;
    
    public:
	Planar();
	void GetUnitVector(unsigned short face_, Vector3 &unit);
	void GetPosition(Vector3 &pos){ pos = position; }
	double GetX(){ return position.axis[0]; }
	double GetY(){ return position.axis[1]; }
	double GetZ(){ return position.axis[2]; }
	double GetWidth(){ return width; }
	double GetLength(){ return length; }
	double GetDepth(){ return depth; }
	bool IsSmall(){ return small; }
	bool IsMedium(){ return med; }
	bool IsLarge(){ return large; }
	void GetLocalCoords(const Vector3&, double&, double&, double&);
	void SetSmall(){ length = 0.60; width = 0.03; depth = 0.03; small = true; med = false; large = false; need_set = true; }
	void SetMedium(){ length = 1.20; width = 0.05; depth = 0.03; small = false; med = true; large = false; need_set = true; }
	void SetLarge(){ length = 2.00; width = 0.05; depth = 0.05; small = false; med = false; large = true; need_set = true; }
	void SetSize(double, double, double);
	void SetPosition(const Vector3&);
	void SetPosition(double, double, double);
	void SetPolarPosition(double, double, double);
	void SetBarRotation(double, double, double);
	void SetUnitVectors(const Vector3&, const Vector3&, const Vector3&);
	bool CheckBounds(unsigned short face_, double x_, double y_, double z_);
	bool PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned short face_, Vector3 &P);
	short FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz);
	short TopBotIntersect(const Vector3 &offset_, const Vector3 &direction, Vector3 &intersect, double &px, double &py, double &pz);
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &P2, short &face1, short &face2);
	double GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned short f1_, unsigned short f2_, Vector3 &intersect1, Vector3 &intersect2);
	std::string DumpVertex();
	std::string DumpDet();
};

class Wall: public Planar{
    private:
	unsigned short num_bars;
	Planar *bars;
	bool init; 
	
    public:
	Wall(){ bars = NULL; init = false; num_bars = 0; }
	~Wall(){ if(init){ delete[] bars; } }
	void Initialize(unsigned short, double, double, double, double);
	unsigned short GetNumBars(){ return num_bars; }
	Planar *GetBar(unsigned short);
	Planar *GetBars(){ return bars; }
	void WallTest();
	std::string DumpDetWall();
};

unsigned short ReadEffFile(const char*, double*, double*);
unsigned short ReadDetFile(const char*, Planar*);
unsigned int TestDetSetup(Planar *bar_array, unsigned short num_bars, unsigned int num_trials);
void DetAng_plan(double, double, unsigned short, double, double, double, bool, bool, double*, double*, double*, double*);
void DetHit_plan(double**, double**, double**, double**, unsigned short, unsigned short*, double*, double, double, unsigned short&, unsigned short&, unsigned short&);
void DetSet_plan_read(const char*, unsigned short, unsigned short*, double*, double*, double*, double*, double*, bool*, bool*, double*, double*, double*, double*, double*, double*, double*, double*);
void det_thick_plan(double, double, double, double&);
void resolution_plan(double, double, double, double, double, double, double, bool, bool, double, double, double, double, unsigned short, unsigned short, double, double, double);
void strag_dE_plan(double, double, double, double, double, double, double, double&, double&, double, double);

#endif
