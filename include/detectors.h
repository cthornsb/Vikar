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
// RegularPolygon
/////////////////////////////////////////////////////////////////////

class RegularPolygon{
  private:
	Line *lines; /// Array of line segments making up the polygon.
	unsigned int nSides; /// The number of sides this polygon contains.
	double chord_length; /// The length of each line segment (m).
	double radius; /// The radius of a circle which bounds the entire polygon (m).
	double sector; /// The angle of each sector of the polygon (rad).
	bool init; /// True if the line array has been initialized.

  public:	
  	/// Default constructor.
	RegularPolygon();
	
	/** Polygon constructor. radius_ is the radius of a circle
	  * which is completely bound within the polygon (in m) and
	  * nSides_ is the number of sides of the polygon.
	  */
	RegularPolygon(const double &radius_, const unsigned int &nSides_){ Initialize(radius_, nSides_); }

	/// Destructor.
	~RegularPolygon(){ if(init){ delete[] lines; } }

	/** Polygon constructor. radius_ is the radius of a circle
	  * which is completely bound within the polygon (in m) and
	  * nSides_ is the number of sides of the polygon.
	  */
	bool Initialize(const double &radius_, const unsigned int &nSides_);
	
	/// Return the number of sides the polygon contains.
	unsigned int GetNsides(){ return nSides; }
	
	/// Return the chord length of each line segment.
	double GetChordLength(){ return chord_length; }
	
	/// Return the radius of a circle which bounds the polygon (m).
	double GetRadius(){ return radius; }
	
	/// Return the angle of each sector (rad);
	double GetSector(){ return sector; }
	
	/** Return true if the point pt_ is contained within the polygon or the
	  * point lies on one of its line segments and return false otherwise.
	  */
	bool IsInside(const Vector3 &pt_);

	/** Return true if the point pt_ is contained within the polygon or the
	  * point lies on one of its line segments and return false otherwise.
	  */	
	bool IsInside(const double &x_, const double &y_);
	
	/// Dump information about the line segments of the polygon.
	void Dump();
};

/////////////////////////////////////////////////////////////////////
// NewVIKARDet
/////////////////////////////////////////////////////////////////////

struct NewVIKARDet{
	float data[9];
	std::string type;
	std::string subtype;
	std::string material;
	unsigned int location;
	
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
	bool need_set; /// True if the global face coordinates need set.

	/** Set the global face coordinates (wrt global origin)
	  * Each vertex is the center coordinate of one of the faces.
	  */
	void _set_face_coords();
	 	
  protected:
    bool use_material; /// True if a material is to be used for energy loss.
    unsigned int material_id; /// The ID of the material to use for energy loss.
	unsigned int front_face; /// The ID of the front face.
	unsigned int back_face; /// The ID of the back face.
	unsigned int location; /// The ID of the detector.
	Vector3 position; /// Center of the detector, cartesian coordinates (in m).
	Vector3 detX, detY, detZ; /// Local detector unit vectors.
	Vector3 GlobalFace[6]; /// The center of the six detector faces in 3d space (in m).
	Matrix3 rotationMatrix; /// The rotation matrix of the detector.
	double length, width, depth; /// Physical size of the detector (in m).
	double theta, phi, psi; /// Rotation of the detector (in radians).
	bool small, med, large; /// The size of the detector if it has type 'vandle'.
	bool use_recoil; /// True if this detector is to be used to detect recoil particles.
	bool use_eject; /// True if this detector is to be used to detect ejectile particles.
	std::string type, subtype; /// The type and subtype of the detector.
	std::string material_name; /// The name of the material to use for energy loss calculations.
    
  public:
    /// Default constructor.
	Planar();
	
	/// Constructor using a NewVIKARDet object.
	Planar(const NewVIKARDet &det_);
	
	virtual ~Planar(){}
	
	/// Return the ID of the material to use for energy loss calculations.
	unsigned int GetMaterial(){ return material_id; }

	/// Return the ID of the "front" face.
	unsigned int GetFrontFace(){ return front_face; }

	/// Return the ID of the "back" face.
	unsigned int GetBackFace(){ return back_face; }

	/// Return the ID of the detector.
	unsigned int GetLoc(){ return location; }

	/// Return the unit vector of one of the faces. Returns the zero vector if the face is undefined
	void GetUnitVector(unsigned int face_, Vector3 &unit);

	/// Assign a vector to the position of the center of the detector in 3d space.
	void GetPosition(Vector3 &pos){ pos = position; }

	/// Get the x-coordinate of the center of the detector (in m).
	double GetX(){ return position.axis[0]; }

	/// Get the y-coordinate of the center of the detector (in m).
	double GetY(){ return position.axis[1]; }

	/// Get the z-coordinate of the center of the detector (in m).
	double GetZ(){ return position.axis[2]; }

	/// Get the width of the detector along the x-axis (in m).
	double GetWidth(){ return width; }

	/// Get the length of the detector along the y-axis (in m).
	double GetLength(){ return length; }

	/// Get the depth of the detector along the z-axis (in m).
	double GetDepth(){ return depth; }

	/// Return the type of the detector.
	std::string GetType(){ return type; }

	/// Return the sub-type of the detector.
	std::string GetSubtype(){ return subtype; }
	
	/// Return the name of the material used for energy loss calculations.
	std::string GetMaterialName(){ return material_name; }

	/// Return the local bar frame coordinates of a global coordinate
	void GetLocalCoords(const Vector3&, double&, double&, double&);
	
	/// Get a vector pointing to a 3d point inside of this geometry
	void GetRandomPointInside(Vector3& output);
	
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

	/// Set the front and rear face ID for this detector.
	void SetFrontFace(unsigned int front_face_);
	
	/// Set the ID of the detector.
	void SetLocation(unsigned int location_){ location = location_; }
	
	void SetMaterialName(std::string name_){ material_name = name_; }

	void SetType(std::string type_){ type = type_; }

	void SetSubtype(std::string subtype_){ subtype = subtype_; }

	void SetSmall(){ length = 0.60; width = 0.03; depth = 0.03; small = true; med = false; large = false; need_set = true; }

	void SetMedium(){ length = 1.20; width = 0.05; depth = 0.03; small = false; med = true; large = false; need_set = true; }

	void SetLarge(){ length = 2.00; width = 0.05; depth = 0.05; small = false; med = false; large = true; need_set = true; }

	/** Set the physical size of the bar
	  * For known bar sizes, it is better to use SetSmall/SetMedium/SetLarge methods
	  * Unknown bar types will not include efficiency data.
	  */
	void SetSize(double length_, double width_, double depth_);
	
	/// Set the position of the center of the bar using a 3d vector (in meters).
	void SetPosition(const Vector3&);
	
	/// Set the cartesian position of the center of the bar in 3d space (in meters).
	void SetPosition(double, double, double);

	/// Set the polar position of the center of the bar in 3d space (meters and radians).
	void SetPolarPosition(double, double, double);
	
	/** Get the local unit vectors using 3d matrix rotation
	  * X and Y are face axes, Z is the axis into or out of the bar.
	  */
	void SetRotation(double, double, double);

	/** Manually set the local bar unit vectors
	  * I would advise against this, use SetRotation instead
	  * Note: This method does not update the bar rotation values
	  *  and should therefore only be used for testing and debugging.
	  */
	void SetUnitVectors(const Vector3 &unitX, const Vector3 &unitY, const Vector3 &unitZ);
	
	void SetRecoil(bool input_=true){ use_recoil = input_; }
	
	void SetEjectile(bool input_=true){ use_eject = input_; }
	
	/** Check if a point (in local coordinates) is within the bounds of the primitive
	  * Return true if the coordinates are within the primitive and false otherwise.
	  */
	virtual bool CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_);

	/** Find if a ray (from the origin) intersects the infinite
	  * plane defined by one of the faces of the VANDLE bar
	  * face 0 is along the +z local axis
	  * face 1 is along the +x local axis
	  * face 2 is along the -z local axis
	  * face 3 is along the -x local axis
	  * face 4 is along the +y local axis
	  * face 5 is along the -y local axis.
	  */
	bool PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, Vector3 &P);
	
	/** Determine if a ray from the origin intersected a face of the bar
	  * and return the face id number. Return -1 if no intersection was found.
	  */
	int FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz);
	
	/** Check for an intersect with the top or bottom faces.
	  * In order to save time, you would normally not want to check
	  * this because the PMTs will block the top and bottom faces.
	  */
	int TopBotIntersect(const Vector3 &offset_, const Vector3 &direction, Vector3 &intersect, double &px, double &py, double &pz);
	
	/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this primitive shape
	  * offset_ is the point where the ray originates wrt the global origin
	  * direction_ is the direction of the ray wrt the global origin
	  * P1 is the first intersection point in global coordinates
	  * P2 is the second intersection point in global coordinates
	  * Return true if the primitive is intersected, and false otherwise.
	  */
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &P2, int &face1, int &face2, double &px, double &py, double &pz);
	
	/** Trace a ray through the detector and calculate the thickness it sees between two faces (f1 and f2)
	  * Return -1 if the ray does not travel through both faces.
	  */
	double GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned int f1_, unsigned int f2_, Vector3 &intersect1, Vector3 &intersect2);
	
	/** Dump raw cartesian face vertex data.
	  * This returns a string containing the vertex coordinates of the
	  * centers of all six faces of the VANDLE bar and its center coordinate.
	  */
	std::string DumpVertex();
	
	/** Dump VIKAR detector format string
	  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)].
	  */
	std::string DumpDet();
};

/////////////////////////////////////////////////////////////////////
// Wall
/////////////////////////////////////////////////////////////////////

class Wall: public Planar{
    private:
	unsigned int num_bars; /// The number of bars in the array wall.
	Planar *bars; /// The array of detectors in the wall.
	bool init; /// True if the wall has been initialized.
	
    public:
    /// Default constructor.
	Wall(){ bars = NULL; init = false; num_bars = 0; }
	
	/// Destructor.
	~Wall(){ if(init){ delete[] bars; } }
	
	/** Initialize the wall object.
	  * Spacing is the physical distance between the edges of bars.
	  * This method leaves a distance of spacing/2 on each side of the wall
	  * so that walls having the same bar spacing will align properly.
	  */
	void Initialize(unsigned int, double, double, double, double);
	
	/// Return the total number of detectors.
	unsigned int GetNumBars(){ return num_bars; }
	
	/** Return a pointer to a bar in this wall
	  * Return null if the bar does not exist.
	  */
	Planar *GetBar(unsigned int);
	
	/// Return a pointer to the array of planar detectors.
	Planar *GetBars(){ return bars; }
	
	/// Test the planar wall class.
	void WallTest();
	
	/** Dump VIKAR detector format string
	  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)].
	  */
	std::string DumpDetWall();
};

/////////////////////////////////////////////////////////////////////
// Elliptical
/////////////////////////////////////////////////////////////////////

class Elliptical : public Planar {
  private:
	double rLong; /// The radius of the long axis of the ellipse (m).
	double rShort; /// The radius of the short axis of the ellipse (m).

  public:
  	/// Default constructor.
	Elliptical() : Planar() {}

	/// Constructor using a NewVIKARDet object.
	Elliptical(const NewVIKARDet &det_) : Planar(det_) { SetRadii(det_.data[6], det_.data[7]); }

	double GetLongAxis() const { return rLong; }
	
	double GetShortAxis() const { return rShort; }

	void SetLongAxis(const double &rLong_){ rLong = rLong_>=0.0?rLong_:-1*rLong_; }

	void SetShortAxis(const double &rShort_){ rShort = rShort_>=0.0?rShort_:-1*rShort_; }

	void SetRadii(const double &rLong_, const double &rShort_);

	/** Check if a point (in local coordinates) is within the bounds of the primitive
	  * Return true if the coordinates are within the primitive and false otherwise.
	  * The Elliptical class checks if the face intersect is within 
	  */
	bool CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_);
};

/////////////////////////////////////////////////////////////////////
// Polygon
/////////////////////////////////////////////////////////////////////

class Polygonal : public Planar {
  private:
	RegularPolygon poly;

  public:
  	/// Default constructor.
	Polygonal() : Planar() {}

	/// Constructor using a NewVIKARDet object.
	Polygonal(const NewVIKARDet &det_) : Planar(det_) { poly.Initialize(det_.data[6], (unsigned int)det_.data[7]); }

	unsigned int GetNsides() { return poly.GetNsides(); }
	
	double GetSectorAngle() { return poly.GetSector(); }
	
	double GetRadius() { return poly.GetRadius(); }
	
	/// Return the chord length of the polygon (m).
	double GetChordLength() { return poly.GetChordLength(); }
	
	/** Check if a point (in local coordinates) is within the bounds of the primitive
	  * Return true if the coordinates are within the primitive and false otherwise.
	  * The Elliptical class checks if the face intersect is within 
	  */
	bool CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_);
};

/////////////////////////////////////////////////////////////////////
// Annular
/////////////////////////////////////////////////////////////////////

class Annular : public Planar {
  private:
	double inRadius; /// The inner radius of the annular detector inside of which no particles will be detected (m).
	double outRadius; /// The outer radius of the annular detector outside of which no particles will be detected (m).
	
  public:
  	/// Default constructor.
	Annular() : Planar() {}

	/// Constructor using a NewVIKARDet object.
	Annular(const NewVIKARDet &det_) : Planar(det_) { SetRadii(det_.data[6], det_.data[7]); }
  
	double GetInnerRadius() const { return inRadius; }
	
	double GetOuterRadius() const { return outRadius; }
	
	void SetInnerRadius(const double &inRadius_){ inRadius = inRadius_>=0.0?inRadius_:-1*inRadius_; }
	
	void SetOuterRadius(const double &outRadius_){ outRadius = outRadius_>=0.0?outRadius_:-1*outRadius_; }

	void SetRadii(const double &inRadius_, const double &outRadius_);
  
	/** Check if a point (in local coordinates) is within the bounds of the primitive
	  * Return true if the coordinates are within the primitive and false otherwise.
	  * The Elliptical class checks if the face intersect is within 
	  */
	bool CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_);
};

/// Read in a detector efficiency file
unsigned int ReadEffFile(const char*, double*, double*);

// Read NewVIKAR detector file and load bars into an array
// Returns the number of detectors loaded from the file
// Assumes the following detector file format for each bar in file
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
int ReadDetFile(const char* fname_, std::vector<Planar*> &detectors);

#endif
