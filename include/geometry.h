/** \file geometry.h
 * \brief Classes which provide support for doing ray-tracing on 3d geometry.
 *
 * Primitive is the main class in this file. All geometry classes (sphere,
 * cylinder, etc.) are derived from the Primitive base class. Primitive
 * simplifies the setup of geometry in 3d space and handles ray tracing
 * operations on said geometry. Primitive may also be used to do rendering
 * of 3d geometry (try make renderer).
 *
 * \author C. R. Thornsberry
 * \date Feb. 26th, 2016
 */
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>

#include "vandmc_core.h"

class NewVIKARdet;

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
// Primitive
/////////////////////////////////////////////////////////////////////

class Primitive{  
  protected:
	bool need_set; /// True if the global face coordinates need set.
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
	bool use_recoil; /// True if this detector is to be used to detect recoil particles.
	bool use_eject; /// True if this detector is to be used to detect ejectile particles.
	bool use_gamma; /// True if this detector is to be used to detect gamma rays.
	bool use_veto; /// True if this detector is to be used to stop particles.
	std::string type, subtype; /// The type and subtype of the detector.
	std::string material_name; /// The name of the material to use for energy loss calculations.

	/** Set the global face coordinates (wrt global origin)
	  * Each vertex is the center coordinate of one of the faces.
	  */
	void _set_face_coords();
    
  public:
    /// Default constructor.
	Primitive();
	
	/// Constructor using a NewVIKARDet object.
	Primitive(NewVIKARdet *det_);
	
	virtual ~Primitive(){}
	
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

	/// Return the local detector frame coordinates of a global coordinate
	void GetLocalCoords(const Vector3&, double&, double&, double&);
	
	/// Get a vector pointing to a 3d point inside of this geometry
	void GetRandomPointInside(Vector3& output);

	bool IsRecoilDet(){ return use_recoil; }

	bool IsEjectileDet(){ return use_eject; }
	
	bool IsGammaDet(){ return use_gamma; }
	
	bool IsVeto(){ return use_veto; }

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

	/// Set the physical size of the detector.
	void SetSize(double length_, double width_, double depth_);
	
	/// Set the position of the center of the detector using a 3d vector (in meters).
	void SetPosition(const Vector3&);
	
	/// Set the cartesian position of the center of the detector in 3d space (in meters).
	void SetPosition(double, double, double);

	/// Set the polar position of the center of the detector in 3d space (meters and radians).
	void SetPolarPosition(double, double, double);
	
	/** Get the local unit vectors using 3d matrix rotation
	  * X and Y are face axes, Z is the axis into or out of the detector.
	  */
	void SetRotation(double, double, double);

	/** Manually set the local detector unit vectors
	  * I would advise against this, use SetRotation instead
	  * Note: This method does not update the detector rotation values
	  *  and should therefore only be used for testing and debugging.
	  */
	void SetUnitVectors(const Vector3 &unitX, const Vector3 &unitY, const Vector3 &unitZ);
	
	void SetRecoil(bool input_=true){ use_recoil = input_; }
	
	void SetEjectile(bool input_=true){ use_eject = input_; }
	
	void SetGamma(bool input_=true){ use_gamma = input_; }
	
	void SetVeto(bool input_=true){ use_gamma = input_; }
	
	/** Check if a point (in local coordinates) is within the bounds of the primitive
	  * Return true if the coordinates are within the primitive and false otherwise.
	  */
	virtual bool CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_);

	/** Find if a ray (from the origin) intersects the infinite
	  * plane defined by one of the faces of the detector
	  * face 0 is along the +z local axis
	  * face 1 is along the +x local axis
	  * face 2 is along the -z local axis
	  * face 3 is along the -x local axis
	  * face 4 is along the +y local axis
	  * face 5 is along the -y local axis.
	  */
	bool PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, double &t);

	/** Find if a ray (from the origin) intersects the infinite cylinder
	  * which bounds this 3d object. The radius of the cylinder is taken
	  * as the "width" of the detector. That is, the size along the x-axis.
	  */
	bool CylinderIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2);

	/** Find if a ray (from the origin) intersects the infinite cone
	  * which bounds this 3d object. The opening angle of the cone is taken
	  * as atan(width/length).
	  */
	bool ConeIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2);

	/** Find if a ray (from the origin) intersects the sphere
	  * which bounds this 3d object. The radius of the bounding
	  * sphere is taken as the "length" of the detector. That is,
	  * the size along the y-axis.
	  */
	bool SphereIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2);
	
	/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this 
	  * primitive shape offset_ is the point where the ray originates wrt the global origin.
	  * direction_ is the direction of the ray wrt the global origin.
	  * t1 is the parameter of the front intersection point given by P1 = offset_ + direction_*t1.
	  * t2 is the parameter of the back intersection point given by P2 = offset_ + direction_*t2.
	  * norm is the normal vector to the surface at point P1.
	  * Return true if the primitive is intersected, and false otherwise.
	  */
	virtual bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2);

	/// Alternate version of IntersectPrimitive which does not return the normal vector.
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, double &t1, double &t2);

	/** Trace a ray through the detector and calculate the thickness it sees between two faces (f1 and f2)
	  * Return -1 if the ray does not travel through both faces.
	  */
	double GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect1, double &t1, double &t2);
	
	/** Dump raw cartesian face vertex data.
	  * This returns a string containing the vertex coordinates of the
	  * centers of all six faces of the VANDLE detector and its center coordinate.
	  */
	std::string DumpVertex();
	
	/** Dump VIKAR detector format string
	  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
	  */
	std::string DumpDet();
};

/////////////////////////////////////////////////////////////////////
// Planar
/////////////////////////////////////////////////////////////////////

class Planar : public Primitive {
  private:
  
  public:
  	/// Default constructor.
	Planar() : Primitive() {}
	
	/// Constructor using a NewVIKARDet object.
	Planar(NewVIKARdet *det_);
	
	void SetSmall(){ SetSize(0.6, 0.03, 0.03); }

	void SetMedium(){ SetSize(1.2, 0.06, 0.03); }

	void SetLarge(){ SetSize(2.0, 0.05, 0.05); }
};

/////////////////////////////////////////////////////////////////////
// Cylindrical
/////////////////////////////////////////////////////////////////////

class Cylindrical : public Primitive { 
  private:
  
  public:
  	/// Default constructor.
	Cylindrical() : Primitive() {}
	
	/// Constructor using a NewVIKARDet object.
	Cylindrical(NewVIKARdet *det_) : Primitive(det_) {}
	
	/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this 
	  * primitive shape offset_ is the point where the ray originates wrt the global origin.
	  * direction_ is the direction of the ray wrt the global origin.
	  * t1 is the parameter of the front intersection point given by P1 = offset_ + direction_*t1.
	  * t2 is the parameter of the back intersection point given by P2 = offset_ + direction_*t2.
	  * norm is the normal vector to the surface at point P1.
	  * Return true if the primitive is intersected, and false otherwise.
	  */
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2);
};

/////////////////////////////////////////////////////////////////////
// Conical
/////////////////////////////////////////////////////////////////////

class Conical : public Primitive { 
  private:
	double alpha; /// The opening angle of the cone.
	double sin_alpha; /// The sine of the opening angle.
	double cos_alpha; /// The cosine of the opening angle.
  
  public:
  	/// Default constructor.
	Conical() : Primitive() {}
	
	/// Constructor using a NewVIKARDet object.
	Conical(NewVIKARdet *det_);
	
	/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this 
	  * primitive shape offset_ is the point where the ray originates wrt the global origin.
	  * direction_ is the direction of the ray wrt the global origin.
	  * t1 is the parameter of the front intersection point given by P1 = offset_ + direction_*t1.
	  * t2 is the parameter of the back intersection point given by P2 = offset_ + direction_*t2.
	  * norm is the normal vector to the surface at point P1.
	  * Return true if the primitive is intersected, and false otherwise.
	  */
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2);
};

/////////////////////////////////////////////////////////////////////
// Spherical
/////////////////////////////////////////////////////////////////////

class Spherical : public Primitive {  
  private:
  
  public:
  	/// Default constructor.
	Spherical() : Primitive() {}
	
	/// Constructor using a NewVIKARDet object.
	Spherical(NewVIKARdet *det_) : Primitive(det_) {}
	
	/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this 
	  * primitive shape offset_ is the point where the ray originates wrt the global origin.
	  * direction_ is the direction of the ray wrt the global origin.
	  * t1 is the parameter of the front intersection point given by P1 = offset_ + direction_*t1.
	  * t2 is the parameter of the back intersection point given by P2 = offset_ + direction_*t2.
	  * norm is the normal vector to the surface at point P1.
	  * Return true if the primitive is intersected, and false otherwise.
	  */
	bool IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2);
};

#endif
