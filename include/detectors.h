// detectors.h
// Cory Thornsberry

#ifndef DETECTORS_H
#define DETECTORS_H

#include <string>

#include "geometry.h"

/////////////////////////////////////////////////////////////////////
// NewVIKARdet
/////////////////////////////////////////////////////////////////////

struct NewVIKARdet{
	float data[9];
	std::string type;
	std::string subtype;
	std::string material;
	unsigned int location;
	
	NewVIKARdet();
	
	NewVIKARdet(std::string input_);
	
	void SetValues(std::string input);
	
	std::string DumpDet();
};

/////////////////////////////////////////////////////////////////////
// Wall
/////////////////////////////////////////////////////////////////////

class Wall: public Primitive{
    private:
	unsigned int num_bars; /// The number of bars in the array wall.
	Primitive *bars; /// The array of detectors in the wall.
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
	Primitive *GetBar(unsigned int);
	
	/// Return a pointer to the array of Primitive detectors.
	Primitive *GetBars(){ return bars; }
	
	/// Test the Primitive wall class.
	void WallTest();
	
	/** Dump VIKAR detector format string
	  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
	  */
	std::string DumpDetWall();
};

/////////////////////////////////////////////////////////////////////
// Elliptical
/////////////////////////////////////////////////////////////////////

class Elliptical : public Primitive {
  private:
	double rLong; /// The radius of the long axis of the ellipse (m).
	double rShort; /// The radius of the short axis of the ellipse (m).

  public:
  	/// Default constructor.
	Elliptical() : Primitive() {}

	/// Constructor using a NewVIKARdet object.
	Elliptical(NewVIKARdet *det_) : Primitive(det_) { SetRadii(det_->data[6], det_->data[7]); }

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

class Polygonal : public Primitive {
  private:
	RegularPolygon poly;

  public:
  	/// Default constructor.
	Polygonal() : Primitive() {}

	/// Constructor using a NewVIKARdet object.
	Polygonal(NewVIKARdet *det_) : Primitive(det_) { poly.Initialize(det_->data[6], (unsigned int)det_->data[7]); }

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

class Annular : public Primitive {
  private:
	double inRadius; /// The inner radius of the annular detector inside of which no particles will be detected (m).
	double outRadius; /// The outer radius of the annular detector outside of which no particles will be detected (m).
	
  public:
  	/// Default constructor.
	Annular() : Primitive() {}

	/// Constructor using a NewVIKARdet object.
	Annular(NewVIKARdet *det_) : Primitive(det_) { SetRadii(det_->data[6], det_->data[7]); }
  
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

/** Read NewVIKAR detector file and load detectors into a vector of pointers.
  * Returns the number of detectors loaded from the file.
  * Assumes the following detector file format for each detector in file
  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
  */
int ReadDetFile(const char* fname_, std::vector<Primitive*> &detectors);

#endif
