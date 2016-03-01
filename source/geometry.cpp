/** \file geometry.cpp
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
#include "geometry.h"
#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// RegularPolygon
/////////////////////////////////////////////////////////////////////

/// Default constructor.
RegularPolygon::RegularPolygon(){
	nSides = 0;
	sector = 0.0;
	radius = 0.0;
	chord_length = 0.0;
	lines = NULL;
	init = false;
}

/** Polygon constructor. radius_ is the radius of a circle
  * which is completely bound within the polygon (in m) and
  * nSides_ is the number of sides of the polygon.
  */
bool RegularPolygon::Initialize(const double &radius_, const unsigned int &nSides_){
	if(init){ return false; }
	
	nSides = nSides_;
	sector = 2.0*pi/nSides;
	radius = radius_/std::cos(sector/2.0);
	chord_length = 2.0*radius*std::sin(sector/2.0);
	
	lines = new Line[nSides];
	
	double theta = -sector/2.0;
	for(unsigned int side = 0; side < nSides; side++){
		lines[side].p1 = Vector3(radius*std::cos(theta), radius*std::sin(theta));
		theta += sector;
		lines[side].p2 = Vector3(radius*std::cos(theta), radius*std::sin(theta));
		lines[side].dir = (lines[side].p2 - lines[side].p1);
	}
	
	return (init = true);
}

/** Return true if the point pt_ is contained within the polygon or the
  * point lies on one of its line segments and return false otherwise.
  */
bool RegularPolygon::IsInside(const Vector3 &pt_){
	return (IsInside(pt_.axis[0], pt_.axis[1]));
}

/** Return true if the point pt_ is contained within the polygon or the
  * point lies on one of its line segments and return false otherwise.
  */
bool RegularPolygon::IsInside(const double &x_, const double &y_){
	Ray trace_ray(x_, y_, x_+1, y_); // Infinite ray along the x-axis.
	Vector3 dummy_vector;
	
	int intersects = 0;
	for(unsigned int side = 0; side < nSides; side++){
		if(lines[side].Intersect(trace_ray, dummy_vector)){ 
			intersects++; 
		}
	}
	
	return (intersects%2!=0?true:false);
}

/** Dump information about the line segments of the polygon.
  *  line# p1x p1y p2x p2y
  */
void RegularPolygon::Dump(){
	for(unsigned int side = 0; side < nSides; side++){
		std::cout << side << "\t" << lines[side].p1.axis[0] << "\t" << lines[side].p1.axis[1] << "\t";
		std::cout << lines[side].p2.axis[0] << "\t" << lines[side].p2.axis[1] << "\n";
	}
}

/////////////////////////////////////////////////////////////////////
// Primitive Class
/////////////////////////////////////////////////////////////////////

// +X is defined as beam-right
// +Y is defined as the vertical axis
// +Z is defined as the beam axis

/// Default constructor.
Primitive::Primitive(){
	material_id = 0;
	front_face = 0;
	back_face = 2;
	location = 0;
	detX = Vector3(1.0, 0.0, 0.0);
	detY = Vector3(0.0, 1.0, 0.0);
	detZ = Vector3(0.0, 0.0, 1.0);
	length = 1.0;
	width = 1.0;
	depth = 1.0;
	need_set = true;
	use_recoil = false;
	use_eject = false;
	use_gamma = false;
	use_veto = false;
	use_material = false;
	type = "unknown";
	subtype = "unknown";
	material_name = "";
}

/// Constructor using a NewVIKARDet object.
Primitive::Primitive(NewVIKARdet *det_){
	material_id = 0;
	front_face = 0;
	back_face = 2;
	location = det_->location;
	detX = Vector3(1.0, 0.0, 0.0);
	detY = Vector3(0.0, 1.0, 0.0);
	detZ = Vector3(0.0, 0.0, 1.0);
	SetPosition(det_->data[0], det_->data[1], det_->data[2]);
	SetRotation(det_->data[3], det_->data[4], det_->data[5]);
	need_set = true;
	SetSize(det_->data[6], det_->data[7], det_->data[8]);
	use_recoil = (det_->type=="dual" || det_->type=="recoil");
	use_eject = (det_->type=="vandle" || det_->type=="dual" || det_->type == "neutron" || det_->type=="eject");
	use_gamma = (det_->type=="vandle" || det_->type == "neutron" || det_->type=="gamma");
	use_veto = (det_->type=="veto");
	use_material = false;
	type = det_->type;
	subtype = det_->subtype;
	material_name = det_->material;	
}

/** Set the global face coordinates (wrt global origin)
  * Each vertex is the center coordinate of one of the faces.
  */
void Primitive::_set_face_coords(){
	// Calculate the center points of each face
	GlobalFace[0] = position-detZ*(depth/2.0); // Front face
	GlobalFace[1] = position+detX*(width/2.0); // Right face
	GlobalFace[2] = position+detZ*(depth/2.0); // Back face
	GlobalFace[3] = position-detX*(width/2.0); // Left face
	GlobalFace[4] = position+detY*(length/2.0); // Top face
	GlobalFace[5] = position-detY*(length/2.0); // Bottom face

	need_set = false;
}

/// Return the unit vector of one of the faces. Returns the zero vector if the face is undefined
void Primitive::GetUnitVector(unsigned int face_, Vector3 &unit){
	if(face_ == 0){ unit = detZ; } // +z local axis face (front)
	else if(face_ == 1){ unit = detX; } // +x local axis (right)
	else if(face_ == 2){ unit = detZ; unit *= -1; } // -z local axis (back)
	else if(face_ == 3){ unit = detX; unit *= -1; } // -x local axis (left)
	else if(face_ == 4){ unit = detY; } // +y local axis (top)
	else if(face_ == 5){ unit = detY; unit *= -1; } // -y local axis (bottom)
	else{ unit = Vector3(0.0, 0.0, 0.0); }
}	

/// Return the local detector frame coordinates of a global coordinate
void Primitive::GetLocalCoords(const Vector3 &world_coords, double &x, double &y, double &z){
	Vector3 temp = (world_coords - position);
	x = temp.Dot(detX);
	y = temp.Dot(detY);
	z = temp.Dot(detZ);
}

/// Get a vector pointing to a 3d point inside of this geometry
void Primitive::GetRandomPointInside(Vector3& output){
	// Get a random point inside the 3d detector
	Vector3 temp = Vector3(frand(-length/2, length/2), frand(-width/2, width/2), frand(-depth/2, depth/2));
	
	// Rotate the random point based on the detector rotation
	rotationMatrix.Transform(temp);
	
	// Get the vector pointing from the origin to the random point
	output = position + temp;
}

/// Set the front and rear face ID for this detector.
void Primitive::SetFrontFace(unsigned int front_face_){ 
	if(front_face_ == 0){ back_face = 2; }
	else if(front_face_ == 1){ back_face = 3; }
	else if(front_face_ == 2){ back_face = 0; }
	else if(front_face_ == 3){ back_face = 1; }
	else if(front_face_ == 4){ back_face = 5; }
	else if(front_face_ == 5){ back_face = 4; }
	else{ return; }
	
	front_face = front_face_;
}

/// Set the physical size of the detector.
void Primitive::SetSize(double length_, double width_, double depth_){
	width = width_; 
	length = length_; 
	depth = depth_;
	need_set = true;
}

/// Set the position of the center of the detector using a 3d vector (in meters).
void Primitive::SetPosition(const Vector3 &pos){
	position = pos;
	need_set = true;
}

/// Set the cartesian position of the center of the detector in 3d space (in meters).
void Primitive::SetPosition(double x, double y, double z){
	position.axis[0] = x; position.axis[1] = y; position.axis[2] = z;
	need_set = true;
}

/// Set the polar position of the center of the detector in 3d space (meters and radians).
void Primitive::SetPolarPosition(double r, double theta, double phi){
	Sphere2Cart(r, theta, phi, position);
	need_set = true;
}

/** Get the local unit vectors using 3d matrix rotation
  * X and Y are face axes, Z is the axis into or out of the detector.
  */
void Primitive::SetRotation(double theta_, double phi_, double psi_){
	theta = theta_; phi = phi_; psi = psi_;
	
	double sin_theta = std::sin(theta); double cos_theta = std::cos(theta);
	double sin_phi = std::sin(phi); double cos_phi = std::cos(phi);
	double sin_psi = std::sin(psi); double cos_psi = std::cos(psi);	
	
	// Pitch-Roll-Yaw convention
	// Rotate by angle theta about the y-axis
	//  angle phi about the z-axis
	//  angle psi about the x-axis
	detX = Vector3(cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta); // Width axis
	detY = Vector3(sin_psi*sin_theta*cos_phi-cos_psi*sin_phi, sin_psi*sin_theta*sin_phi+cos_psi*cos_phi, cos_theta*sin_psi); // Length axis
	detZ = Vector3(cos_psi*sin_theta*cos_phi+sin_psi*sin_phi, cos_psi*sin_theta*sin_phi-sin_psi*cos_phi, cos_theta*cos_psi); // Depth axis
	
	// Normalize the unit vectors
	detX.Normalize(); detY.Normalize(); detZ.Normalize();
	need_set = true;
	
	// Set the rotation matrix
	rotationMatrix.SetUnitX(detX);
	rotationMatrix.SetUnitY(detY);
	rotationMatrix.SetUnitZ(detZ);
}

/** Manually set the local detector unit vectors
  * I would advise against this, use SetRotation instead
  * Note: This method does not update the detector rotation values
  *  and should therefore only be used for testing and debugging.
  */
void Primitive::SetUnitVectors(const Vector3 &unitX, const Vector3 &unitY, const Vector3 &unitZ){
	detX = unitX; detY = unitY; detZ = unitZ;
	
	detX.Normalize(); detY.Normalize(); detZ.Normalize();
	need_set = true;
	
	// Set the rotation matrix
	rotationMatrix.SetUnitX(detX);
	rotationMatrix.SetUnitY(detY);
	rotationMatrix.SetUnitZ(detZ);
}

/** Check if a point (in local coordinates) is within the bounds of the primitive
  * Return true if the coordinates are within the primitive and false otherwise.
  */
bool Primitive::CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_){
	if(face_ == 0 || face_ == 2){ // Front face (upstream) or back face (downstream)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 1 || face_ == 3){ // Right face (beam-right) or left face (beam-left)
		if((z_ >= -depth/2.0 && z_ <= depth/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 4 || face_ == 5){ // Top face (+y) or bottom face (-y)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (z_ >= -depth/2.0 && z_ <= depth/2.0)){ return true; }
	}
	
	return false;	
}

/** Find if a ray (from the origin) intersects the infinite
  * plane defined by one of the faces of the VANDLE detector
  * face 0 is along the +z local axis
  * face 1 is along the +x local axis
  * face 2 is along the -z local axis
  * face 3 is along the -x local axis
  * face 4 is along the +y local axis
  * face 5 is along the -y local axis.
  */
bool Primitive::PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, double &t){
	if(need_set){ _set_face_coords(); }
	
	Vector3 unit;
	GetUnitVector(face_, unit);
	
	t = (GlobalFace[face_].Dot(unit)-offset_.Dot(unit))/(direction_.Dot(unit));
	
	if(t >= 0){ return true; } // The ray intersects the plane
	return false; // The plane is behind the ray, the ray will never intersect it
}

/** Find if a ray (from the origin) intersects the infinite cylinder
  * which bounds this 3d object. The radius of the cylinder is taken
  * as the "width" of the detector. That is, the size along the x-axis.
  */
bool Primitive::CylinderIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	// position is the position of the center of the detector.
	// The direction of the cylinder is given by the local unit vector detY.
	Vector3 dP = offset_ - position;
	Vector3 A1 = direction_ - detY*detY.Dot(direction_);
	Vector3 C1 = dP - detY*detY.Dot(dP);
	
	double A2 = A1.Square();
	double B2 = 2.0 * C1.Dot(A1);
	double C2 = C1.Square() - width*width/4.0;
	
	t1 = (-B2 + std::sqrt(B2*B2-4.0*A2*C2))/(2.0*A2);
	t2 = (-B2 - std::sqrt(B2*B2-4.0*A2*C2))/(2.0*A2);
	
	if(!std::isnan(t1) && !std::isnan(t2)){ // Check that a solution exists.
		if(((t1 > 0 && t2 > 0) && t1 > t2) || t2 > 0){ // Swap the two points since t2 is closer to the ray origin.
			double temp = t1;
			t1 = t2;
			t2 = temp;
		}
		return (t1 >= 0 || t2 >= 0);
	}
	
	return false;
}

/** Find if a ray (from the origin) intersects the infinite cone
  * which bounds this 3d object. The opening angle of the cone is taken
  * as atan(width/length).
  */
bool Primitive::ConeIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	Vector3 dP = offset_ - position;
	Vector3 A1 = direction_ - detY*detY.Dot(direction_);
	double B1 = direction_.Dot(detY);
	Vector3 C1 = dP - detY*detY.Dot(dP);
	double D1 = dP.Dot(detY);
	
	// The squared cosine and sine of the opening angle of the cone.
	double cos2_alpha = length*length / (length*length + width*width/4.0);
	double sin2_alpha = width*width / (4.0*length*length + width*width);
	
	double A2 = cos2_alpha*A1.Square() - sin2_alpha*B1*B1;
	double B2 = 2*cos2_alpha*A1.Dot(C1) - 2*sin2_alpha*B1*D1;
	double C2 = cos2_alpha*C1.Square() - sin2_alpha*D1*D1;
	
	t1 = (-B2 + std::sqrt(B2*B2-4.0*A2*C2))/(2.0*A2);
	t2 = (-B2 - std::sqrt(B2*B2-4.0*A2*C2))/(2.0*A2);
	
	if(!std::isnan(t1) && !std::isnan(t2)){ // Check that a solution exists.
		if(((t1 > 0 && t2 > 0) && t1 > t2) || t2 > 0){ // Swap the two points since t2 is closer to the ray origin.
			double temp = t1;
			t1 = t2;
			t2 = temp;
		}
		return (t1 >= 0 || t2 >= 0);
	}
	
	return false;
}

/** Find if a ray (from the origin) intersects the sphere
  * which bounds this 3d object. The radius of the bounding
  * sphere is taken as the "length" of the detector. That is,
  * the size along the y-axis.
  */
bool Primitive::SphereIntersect(const Vector3 &offset_, const Vector3 &direction_, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	Vector3 R = offset_ - position; // Vector pointing from the center of the detector to the start position of the ray.
	double A = direction_.Square();
	double B = 2.0 * R.Dot(direction_); 
	double C = R.Square() - length*length/4.0;
	
	t1 = (-B + std::sqrt(B*B-4.0*A*C))/(2.0*A);
	t2 = (-B - std::sqrt(B*B-4.0*A*C))/(2.0*A);
	
	if(!std::isnan(t1) && !std::isnan(t2)){ // Check that a solution exists.
		if(((t1 > 0 && t2 > 0) && t1 > t2) || t2 > 0){ // Swap the two points since t2 is closer to the ray origin.
			double temp = t1;
			t1 = t2;
			t2 = temp;
		}
		return (t1 >= 0 || t2 >= 0);
	}
	
	return false;
}

/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this prism. 
  * offset_ is the point where the ray originates wrt the global origin.
  * direction_ is the direction of the ray wrt the global origin.
  * P1 is the first intersection point in global coordinates.
  * P2 is the second intersection point in global coordinates.
  * norm is the normal vector to the surface at point P1.
  * Return true if the prism is intersected, and false otherwise.
  */
bool Primitive::IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	int face_count = 0;
	int face1 = -1;
	double temp_t;
	double px, py, pz;
	for(unsigned int i = 0; i < 6; i++){
		if(PlaneIntersect(offset_, direction_, i, temp_t)){ // The ray intersects the plane containing this face
			// Transform the intersection point into local coordinates and check if they're within the bounds
			GetLocalCoords((offset_ + direction_*temp_t), px, py, pz);
			if(CheckBounds(i, px, py, pz)){ // The face was struck
				face_count++;
				if(face_count == 1){ 
					t1 = temp_t; 
					face1 = i;
				}
				else if(face_count == 2){ // Should not have more than 2 face intersects. Do this just to save a little time.
					if(t1 > temp_t){ // Swap the two points since t2 is closer to the ray origin.
						t2 = t1;
						t1 = temp_t;
						face1 = i;
					}
					else{ t2 = temp_t; }
					break; 
				} 
			}
		} // if(PlaneIntersect(offset_, direction_, i, ray))
	} // for(unsigned int i = 0; i < 6; i++)

	if(face_count == 0){ return false; }
	
	P1 = offset_ + direction_*t1;
	
	// Get the surface normal at point P1.
	GetUnitVector(face1, norm);
	
	return true; 
}

/// Alternate version of IntersectPrimitive which does not return the normal vector.
bool Primitive::IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, double &t1, double &t2){
	Vector3 dummy;
	return IntersectPrimitive(offset_, direction_, P1, dummy, t1, t2);
}
	
/** Trace a ray through the object and calculate the thickness it sees.
  * Return -1 if the ray does not travel through both faces.
  */
double Primitive::GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect1, double &t1, double &t2){
	Vector3 dummy;

	// Check that the ray travels through the object.
	if(!IntersectPrimitive(offset_, direction_, intersect1, dummy, t1, t2)){ return -1; }
	
	return (t2-t1)*direction_.Length();
}

/** Dump raw cartesian face vertex data.
  * This returns a string containing the vertex coordinates of the
  * centers of all six faces of the VANDLE detector and its center coordinate.
  */
std::string Primitive::DumpVertex(){
	if(need_set){ _set_face_coords(); }
	std::stringstream stream;
	for(unsigned int i = 0; i < 6; i++){
		stream << GlobalFace[i].axis[0] << "\t" << GlobalFace[i].axis[1] << "\t" << GlobalFace[i].axis[2] << "\n";
	}
	stream << position.axis[0] << "\t" << position.axis[1] << "\t" << position.axis[2];
	return stream.str();
}

/** Dump VIKAR detector format string
  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
  */
std::string Primitive::DumpDet(){
	if(need_set){ _set_face_coords(); }
	std::stringstream stream;
	stream << position.axis[0] << "\t" << position.axis[1] << "\t" << position.axis[2];
	stream << "\t" << theta << "\t" << phi << "\t" << psi;
	stream << "\t" << type << "\t" << subtype;
	if(type != "vandle"){
		stream << "\t" << length << "\t" << width << "\t" << depth;
		stream << "\t" << material_name;
	}
	return stream.str();
}

/////////////////////////////////////////////////////////////////////
// Planar
/////////////////////////////////////////////////////////////////////

/// Constructor using a NewVIKARDet object.
Planar::Planar(NewVIKARdet *det_) : Primitive(det_) {
	if(det_->type == "vandle"){
		if(det_->subtype == "small"){ SetSmall(); }
		else if(det_->subtype == "medium"){ SetMedium(); }
		else if(det_->subtype == "large"){ SetLarge(); }
		else{ 
			std::cout << " Warning! Unknown VANDLE subtype = " << det_->subtype << "!\n";
			SetSize(det_->data[6], det_->data[7], det_->data[8]);			
		}
	}
	else{ SetSize(det_->data[6], det_->data[7], det_->data[8]); }
}

/////////////////////////////////////////////////////////////////////
// Cylindrical
/////////////////////////////////////////////////////////////////////

/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this cylinder.
  * offset_ is the point where the ray originates wrt the global origin.
  * direction_ is the direction of the ray wrt the global origin.
  * P1 is the first intersection point in global coordinates.
  * P2 is the second intersection point in global coordinates.
  * norm is the normal vector to the surface at point P1.
  * Return true if the cylinder is intersected, and false otherwise.
  */
bool Cylindrical::IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	// Check if the ray intersects the infinite cylinder with r = width/2.0.
	if(!CylinderIntersect(offset_, direction_, t1, t2)){ return false; }
	
	P1 = offset_ + direction_*t1;
	
	Vector3 r1 = P1 - position;
	int face1 = 0;
	
	// Check if the intersection points are within the cylinder endcaps.
	double z1 = std::sqrt(r1.Square() - width*width/4.0);
	
	if(r1.CosAngle(detY) < 0.0){ z1 *= -1; }
	
	bool check1 = fabs(z1) <= length/2.0;
	
	if(!check1){ // Point P1 is outside of the bounds of the cylinder.
		double temp_t;
		double px, py, pz;
		for(int i = 4; i < 6; i++){ // Check for intersections with the endcaps. Faces 4 and 5 represent the endcaps of the cylinder.
			if(PlaneIntersect(offset_, direction_, i, temp_t)){
				// Transform the intersection point into local coordinates and check if they're within the bounds
				GetLocalCoords((offset_ + direction_*temp_t), px, py, pz);
				if((px*px + pz*pz) <= width*width/4.0){ // The face was struck
					P1 = offset_ + direction_*temp_t; 
					t1 = temp_t;
					face1 = i;
					break;
				}
			}
		}
	}
	
	if(face1 == 0){ // Calculate the normal to the curved surface at point P1.
		if(!check1){ return false; } // Found no intersection with the endcaps.
		norm = P1 - position - detY*z1;
		norm.Normalize();
	}
	else{ // Get the normal of one of the endcaps.
		GetUnitVector(face1, norm); 
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////
// Conical
/////////////////////////////////////////////////////////////////////

/// Constructor using a NewVIKARDet object.
Conical::Conical(NewVIKARdet *det_) : Primitive(det_) {
	alpha = std::atan2(width/2.0, length); 
	sin_alpha = std::sin(alpha);
	cos_alpha = std::cos(alpha);
}

/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this cone.
  * offset_ is the point where the ray originates wrt the global origin.
  * direction_ is the direction of the ray wrt the global origin.
  * P1 is the first intersection point in global coordinates.
  * P2 is the second intersection point in global coordinates.
  * norm is the normal vector to the surface at point P1.
  * Return true if the cylinder is intersected, and false otherwise.
  */
bool Conical::IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	// Check if the ray intersects the infinite cylinder with r = width/2.0.
	if(!ConeIntersect(offset_, direction_, t1, t2)){ return false; }
	
	P1 = offset_ + direction_*t1;
	
	Vector3 r1 = P1 - position;
	if(r1.CosAngle(detY) > 0.0){ return false; }
	int face1 = 0;

	// Check if the intersection points are within the cone.
	double z1 = r1.Length()*cos_alpha;
	
	bool check1 = fabs(z1) <= length;
	
	if(!check1){ // Point P1 is outside of the bounds of the cylinder.
		double temp_t;
		double px, py, pz;
		if(PlaneIntersect(offset_, direction_, 5, temp_t)){
			// Transform the intersection point into local coordinates and check if they're within the bounds
			GetLocalCoords((offset_ + direction_*temp_t), px, py, pz);
			if((px*px + pz*pz) <= width*width/4.0){ // The face was struck
				t1 = temp_t;
				face1 = 5;
			}
		}
	}

	if(face1 == 0){ // Calculate the normal to the curved surface at point P1.
		if(!check1){ return false; } // Found no intersection with the endcaps.
		Vector3 rprime = P1 - position + detY*z1;
		rprime.Normalize();
		norm = rprime*cos_alpha + detY*sin_alpha;
	}
	else{ // Get the normal of one of the endcaps.
		GetUnitVector(face1, norm); 
	}
		
	return true;
}

/////////////////////////////////////////////////////////////////////
// Spherical
/////////////////////////////////////////////////////////////////////

/** Calculate the intersection of a ray of the form (offset_ + t * direction_) with this sphere.
  * offset_ is the point where the ray originates wrt the global origin.
  * direction_ is the direction of the ray wrt the global origin.
  * P1 is the first intersection point in global coordinates.
  * P2 is the second intersection point in global coordinates.
  * norm is the normal vector to the surface at point P1.
  * Return true if the sphere is intersected, and false otherwise.
  */
bool Spherical::IntersectPrimitive(const Vector3& offset_, const Vector3& direction_, Vector3 &P1, Vector3 &norm, double &t1, double &t2){
	if(need_set){ _set_face_coords(); }
	
	// Check if the ray intersects the sphere with r = length/2.0.
	if(!SphereIntersect(offset_, direction_, t1, t2)){ return false; }
	
	// Calculate the close intersection point.
	P1 = offset_ + direction_*t1;	
	
	// Calculate the normal to the surface at point P1.
	norm = P1 - position;
	norm.Normalize();
	
	return true;
}
