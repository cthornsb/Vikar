// Det_Planar.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:20 2014

#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////////

const Vector3 zero_vector = Vector3(0.0, 0.0, 0.0);

/////////////////////////////////////////////////////////////////////
// Planar Class
//  Used to generate and control individual VANDLE bars
/////////////////////////////////////////////////////////////////////

// +X is defined as beam-right
// +Y is defined as the vertical axis
// +Z is defined as the beam axis

// Default constructor
Planar::Planar(){
	barX = Vector3(1.0, 0.0, 0.0);
	barY = Vector3(0.0, 1.0, 0.0);
	barZ = Vector3(0.0, 0.0, 1.0);
	length = 1.0;
	width = 1.0;
	depth = 1.0;
	small = false;
	med = false;
	large = false;
	need_set = true;
}

// Set the global face coordinates (wrt global origin)
// Each vertex is the center coordinate of one of the faces 
void Planar::_set_face_coords(){
	// Calculate the center points of each face
	GlobalFace[0] = position-barZ*(depth/2.0); // Front face
	GlobalFace[1] = position+barX*(width/2.0); // Right face
	GlobalFace[2] = position+barZ*(depth/2.0); // Back face
	GlobalFace[3] = position-barX*(width/2.0); // Left face
	GlobalFace[4] = position+barY*(length/2.0); // Top face
	GlobalFace[5] = position-barY*(length/2.0); // Bottom face
	
	// Calculate the vertices for all six faces
	// This array will contain duplicate data for simplified calculations later
	// Vertices are indexed from top left clockwise about the normal of the face
	/*GlobalVertex[0][0] = GlobalFace[0] - barX*(width/2.0) + barY*(length/2.0); // -x +y
	GlobalVertex[0][1] = GlobalFace[0] + barX*(width/2.0) + barY*(length/2.0); // +x +y
	GlobalVertex[0][2] = GlobalFace[0] + barX*(width/2.0) - barY*(length/2.0); // +x -y
	GlobalVertex[0][3] = GlobalFace[0] - barX*(width/2.0) - barY*(length/2.0); // -x -y

	GlobalVertex[2][0] = GlobalFace[2] + barX*(width/2.0) + barY*(length/2.0); // +x +y
	GlobalVertex[2][1] = GlobalFace[2] + barX*(width/2.0) - barY*(length/2.0); // +x -y
	GlobalVertex[2][2] = GlobalFace[2] - barX*(width/2.0) - barY*(length/2.0); // -x -y
	GlobalVertex[2][3] = GlobalFace[2] - barX*(width/2.0) + barY*(length/2.0); // -x +y
	
	GlobalVertex[1] = {GlobalVertex[0][1], GlobalVertex[1][0], GlobalVertex[1][3], GlobalVertex[0][2]}; // +y -z | +y +z | -y +z | -y -z
	GlobalVertex[3] = {GlobalVertex[1][1], GlobalVertex[0][0], GlobalVertex[0][3], GlobalVertex[1][2]}; // +y +z | +y -z | -y -z | -y +z
	GlobalVertex[4] = {GlobalVertex[1][1], GlobalVertex[1][0], GlobalVertex[0][1], GlobalVertex[0][0]}; // -x +z | +x +z | +x -z | -x -z
	GlobalVertex[5] = {GlobalVertex[1][3], GlobalVertex[1][2], GlobalVertex[0][3], GlobalVertex[0][2]}; // -x +z | +x +z | +x -z | -x -z*/

	need_set = false;
}

// Return the unit vector of one of the faces. Returns the zero vector if the face is undefined
void Planar::GetUnitVector(unsigned short face_, Vector3 &unit){
	if(face_ == 0){ unit = barZ; } // +z local axis face (front)
	else if(face_ == 1){ unit = barX; } // +x local axis (right)
	else if(face_ == 2){ unit = barZ; unit *= -1; } // -z local axis (back)
	else if(face_ == 3){ unit = barX; unit *= -1; } // -x local axis (left)
	else if(face_ == 4){ unit = barY; } // +y local axis (top)
	else if(face_ == 5){ unit = barY; unit *= -1; } // -y local axis (bottom)
	else{ unit = Vector3(0.0, 0.0, 0.0); }
}	

// Return the local bar frame coordinates of a global coordinate
void Planar::GetLocalCoords(const Vector3 &world_coords, double &x, double &y, double &z){
	Vector3 temp = (world_coords - position);
	x = temp.Dot(barX);
	y = temp.Dot(barY);
	z = temp.Dot(barZ);
}

// Set the physical size of the bar
// For known bar sizes, it is better to use SetSmall/SetMedium/SetLarge methods
// Unknown bar types will not include efficiency data
void Planar::SetSize(double length_, double width_, double depth_){
	if(length_ == 0.6 && width_ == 0.03 && depth_ == 0.03){ SetSmall(); }
	else if(length_ == 1.2 && width_ == 0.05 && depth_ == 0.03){ SetMedium(); }
	else if(length_ == 2.0 && width_ == 0.05 && depth_ == 0.05){ SetLarge(); }
	else{
		width = width_; length = length_; depth = depth_; 
		small = false; med = false; large = false; need_set = true;
	}
}

// Set the position of the center of the bar using a 3d vector (in meters)
void Planar::SetPosition(const Vector3 &pos){
	position = pos;
	need_set = true;
}

// Set the cartesian position of the center of the bar in 3d space (in meters)
void Planar::SetPosition(double x, double y, double z){
	position.axis[0] = x; position.axis[1] = y; position.axis[2] = z;
	need_set = true;
}

// Set the polar position of the center of the bar in 3d space (meters and radians)
void Planar::SetPolarPosition(double r, double theta, double phi){
	Sphere2Cart(r, theta, phi, position);
	need_set = true;
}

// Get the local unit vectors using 3d matrix rotation
// X and Y are face axes, Z is the axis into or out of the bar
void Planar::SetBarRotation(double theta_, double phi_, double psi_){
	theta = theta_; phi = phi_; psi = psi_;

	/*double sin_theta = std::sin(theta); double cos_theta = std::cos(theta);
	double sin_phi = std::sin(phi); double cos_phi = std::cos(phi);
	double sin_psi = std::sin(psi); double cos_psi = std::cos(psi);
	
	barX = Vector3(cos_psi*cos_phi-cos_theta*sin_phi*sin_psi, cos_psi*sin_phi+cos_theta*cos_phi*sin_psi, sin_psi*sin_theta); // Width axis
	barY = Vector3(-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi, -sin_psi*sin_phi+cos_theta*cos_phi*cos_psi, cos_psi*sin_theta); // Length axis
	barZ = Vector3(sin_theta*sin_phi, -sin_theta*cos_phi, cos_theta);  // Depth axis*/
	
	double sin_theta = std::sin(theta); double cos_theta = std::cos(theta);
	double sin_phi = std::sin(phi); double cos_phi = std::cos(phi);
	double sin_psi = std::sin(psi); double cos_psi = std::cos(psi);
	
	// Pitch-Roll-Yaw convention
	// Rotate by angle theta about the y-axis
	//  angle phi about the z-axis
	//  angle psi about the x-axis
	barX = Vector3(cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta); // Width axis
	barY = Vector3(sin_psi*sin_theta*cos_phi-cos_psi*sin_phi, sin_psi*sin_theta*sin_phi+cos_psi*cos_phi, cos_theta*sin_psi); // Length axis
	barZ = Vector3(cos_psi*sin_theta*cos_phi+sin_psi*sin_phi, cos_psi*sin_theta*sin_phi-sin_psi*cos_phi, cos_theta*cos_psi); // Depth axis
	
	// Normalize the unit vectors
	barX.Normalize(); barY.Normalize(); barZ.Normalize();
	need_set = true;    		
}

// Manually set the local bar unit vectors
// I would advise against this, use SetBarRotation instead
// Note: This method does not update the bar rotation values
//  and should therefore only be used for testing and debugging
void Planar::SetUnitVectors(const Vector3 &unitX, const Vector3 &unitY, const Vector3 &unitZ){
	barX = unitX; barY = unitY; barZ = unitZ;
	barX.Normalize(); barY.Normalize(); barZ.Normalize();
	need_set = true;
}

// Check if a point (in local coordinates) is within the bounds of the primitive
// Return true if the coordinates are within the primitive and false otherwise
bool Planar::CheckBounds(unsigned short face_, double x_, double y_, double z_){
	if(face_ == 0){ // Front face (upstream)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 1){ // Right face (beam-right)
		if((z_ >= -depth/2.0 && z_ <= depth/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 2){ // Back face (downstream)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 3){ // Left face (beam-left)
		if((z_ >= -depth/2.0 && z_ <= depth/2.0) && (y_ >= -length/2.0 && y_ <= length/2.0)){ return true; }
	}
	else if(face_ == 4){ // Top face (+y)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (z_ >= -depth/2.0 && z_ <= depth/2.0)){ return true; }
	}
	else if(face_ == 5){ // Bottom face (-y)
		if((x_ >= -width/2.0 && x_ <= width/2.0) && (z_ >= -depth/2.0 && z_ <= depth/2.0)){ return true; }
	}
	return false;	
}

// Find if a ray (from the origin) intersects the infinite
// plane defined by one of the faces of the VANDLE bar
// face 0 is along the +z local axis
// face 1 is along the +x local axis
// face 2 is along the -z local axis
// face 3 is along the -x local axis
// face 4 is along the +y local axis
// face 5 is along the -y local axis
bool Planar::PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned short face_, Vector3 &P){
	if(need_set){ _set_face_coords(); }
	
	Vector3 unit;
	GetUnitVector(face_, unit);
	
	// The ray vector has the parametric form ray = offset_ + t*direction_
	// First find the intersection points between the ray and a plane containing the face polygon
	double t = (GlobalFace[face_].Dot(unit)-offset_.Dot(unit))/(direction_.Dot(unit));
	P = offset_ + direction_*t; // The intersection point on the plane
	
	if(t >= 0){ return true; } // The ray intersects the plane
	return false; // The plane is behind the ray, the ray will never intersect it
}

// Determine if a ray from the origin intersected a face of the bar
// and return the face id number. Return -1 if no intersection was found
short Planar::FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz){
	if(need_set){ _set_face_coords(); }
	
	// Check the Front face (Normal to local +z axis)
	if(PlaneIntersect(offset_, direction_, 0, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (py >= -length/2.0 && py <= length/2.0)){ 
			return 0; // Ray intersects the front face				
		}
	}

	// Check the Right face (Normal to local +x axis)
	if(PlaneIntersect(offset_, direction_, 1, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((pz >= -depth/2.0 && pz <= depth/2.0) && (py >= -length/2.0 && py <= length/2.0)){
			return 1; // Ray intersects the right face
		}
	}

	// Check the Back face (Normal to local -z axis)
	if(PlaneIntersect(offset_, direction_, 2, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (py >= -length/2.0 && py <= length/2.0)){ 
			return 2; // Ray intersects the back face     				
		}
	}

	// Check the Left face (Normal to local -x axis)
	if(PlaneIntersect(offset_, direction_, 3, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((pz >= -depth/2.0 && pz <= depth/2.0) && (py >= -length/2.0 && py <= length/2.0)){
			return 3; // Ray intersects the left face
		}
	}
	
	return -1; // No intersection
}

// Check for an intersect with the top or bottom faces.
// In order to save time, you would normally not want to check
// this because the PMTs will block the top and bottom faces
short Planar::TopBotIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz){
	if(need_set){ _set_face_coords(); }
	
	// Check the Top face (Normal to local +y axis)
	if(PlaneIntersect(offset_, direction_, 4, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (pz >= -depth/2.0 && pz <= depth/2.0)){ 
			return 4; // Ray intersects the top face				
		}
	}
	
	// Check the Bottom face (Normal to local -y axis)
	if(PlaneIntersect(offset_, direction_, 5, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (pz >= -depth/2.0 && pz <= depth/2.0)){ 
			return 5; // Ray intersects the bottom face				
		}
	}
	
	return -1;
}

// Calculate the intersection of a ray of the form (offset_ + t * direction_) with this primitive shape
// offset_ is the point where the ray originates wrt the global origin
// direction_ is the direction of the ray wrt the global origin
// P1 is the first intersection point in global coordinates
// P2 is the second intersection point in global coordinates
// Return true if the primitive is intersected, and false otherwise
bool Planar::IntersectPrimitive(const Vector3 &offset_, const Vector3 &direction_, Vector3 &P1, Vector3 &P2, short &face1, short &face2){
	if(need_set){ _set_face_coords(); }
	
	short face_count = 0;
	double px, py, pz;
	Vector3 ray, unit;
	for(unsigned short i = 0; i < 6; i++){
		GetUnitVector(i, unit);
		if(PlaneIntersect(offset_, direction_, i, ray)){ // The ray intersects the plane containing this face
			// Transform the intersection point into local coordinates and check if they're within the bounds
			GetLocalCoords(ray, px, py, pz);
			if(CheckBounds(i, px, py, pz)){ // The face was struck
				if(face_count == 0){ 
					P1 = ray; 
					face1 = i; 
				}
				else if(face_count == 1){ 
					P2 = ray; 
					face2 = i; 
					//break;
				}
				face_count++;
			}
		} // if(PlaneIntersect(offset_, direction_, i, ray))
	} // for(unsigned short i = 0; i < 6; i++)
	
	if(face_count == 2){ return true; }
	else if(face_count == 0){ return false; }
	else{ std::cout << " Unexpected number of face intersects: " << face_count << std::endl; }
	return false; 
}

// Trace a ray through the detector and calculate the thickness it sees between two faces (f1 and f2)
// Return -1 if the ray does not travel through both faces
double Planar::GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned short f1_, unsigned short f2_, Vector3 &intersect1, Vector3 &intersect2){
	if(f1_ > 5 || f2_ > 5){ return -1; } // Invalid face number
	if(need_set){ _set_face_coords(); }
	double px, py, pz;
	
	// Check face 1
	if(PlaneIntersect(offset_, direction_, f1_, intersect1)){ 
		GetLocalCoords(intersect1, px, py, pz);
		if(!CheckBounds(f1_, px, py, pz)){ return -1; }
	}
	else{ return -1; }
	
	// Check face 2
	if(PlaneIntersect(offset_, direction_, f2_, intersect2)){
		GetLocalCoords(intersect2, px, py, pz);
		if(!CheckBounds(f2_, px, py, pz)){ return -1; }
	}
	else{ return -1; }
	return Dist3d(intersect1, intersect2);
}

// Dump raw cartesian face vertex data.
// This returns a string containing the vertex coordinates of the
// centers of all six faces of the VANDLE bar and its center coordinate
std::string Planar::DumpVertex(){
	if(need_set){ _set_face_coords(); }
	std::string output = "";
	for(unsigned short i = 0; i < 6; i++){
		output += to_str(GlobalFace[i].axis[0]) + "\t" + to_str(GlobalFace[i].axis[1]) + "\t" + to_str(GlobalFace[i].axis[2]) + "\n";
	}
	output += to_str(position.axis[0]) + "\t" + to_str(position.axis[1]) + "\t" + to_str(position.axis[2]);
	return output;
}

// Dump VIKAR detector format string
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
std::string Planar::DumpDet(){
	if(need_set){ _set_face_coords(); }
	std::string output = to_str(position.axis[0]) + "\t" + to_str(position.axis[1]) + "\t" + to_str(position.axis[2]);
	output += "\t" + to_str(theta) + "\t" + to_str(phi) + "\t" + to_str(psi);
	if(large){ output += "\tlarge\n"; }
	else if(med){ output += "\tmedium\n"; }
	else if(small){ output += "\tsmall\n"; }
	else{ output += "\tunknown\t" + to_str(length) + "\t" + to_str(width) + "\t" + to_str(depth) + "\n"; }
	return output;
}

/////////////////////////////////////////////////////////////////////
// Wall class
//  Used to generate flat walls of bars
/////////////////////////////////////////////////////////////////////

// Initialize the wall object.
// Spacing is the physical distance between the edges of bars.
// This method leaves a distance of spacing/2 on each side of the wall
// so that walls having the same bar spacing will align properly
void Wall::Initialize(unsigned short num_bars_, double spacing, double bar_length, double bar_width, double bar_depth){
	// Bar spacing for Large bar walls = 0.01350 m (4 bars per 10" = 0.254 m)
	// Bar spacing for Small bar ribs = 0.01488 m (7 bars per 36 deg @ 0.5 m)
	if(init){ 
		std::cout << " Warning: Wall already initialized\n";
		return; 
	}  
	else if(num_bars_ == 0){
		std::cout << " Warning: Cannot initialize zero bars\n";
		return; 
	}
		
	num_bars = num_bars_;
	bars = new Planar[num_bars];

	length = bar_length;
	width = num_bars*(spacing + bar_width);  
	depth = bar_depth;
	for(unsigned short i = 0; i < num_bars; i++){
		bars[i].SetSize(bar_length, bar_width, bar_depth);
		bars[i].SetBarRotation(theta, phi, psi);
		bars[i].SetPosition(position + barX*((-width/2.0+spacing/2.0+bar_width/2.0)+i*(bar_width+spacing)));
	}

	init = true;
}

// Return a pointer to a bar in this wall
// Return null if the bar does not exist
Planar* Wall::GetBar(unsigned short bar){ 
	if(init && bar < num_bars){ return &bars[bar]; }
	return NULL;
}

// Test the planar wall class
void Wall::WallTest(){
	if(!init){ std::cout << " Initialized: No\n"; return; }
	else{ std::cout << " Initialized: Yes\n"; }
	
	std::cout << " Length: " << length << std::endl;
	std::cout << " Width: " << width << std::endl;
	std::cout << " Depth: " << depth << std::endl;
}

// Dump VIKAR detector format string
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
std::string Wall::DumpDetWall(){
	std::string output = "";
	if(!init){ 
		std::cout << " Warning: Wall is un-initialized\n";
		return output; 
	}
	
	// Get detector strings for all detectors in wall
	for(unsigned short i = 0; i < num_bars; i++){
		output += bars[i].DumpDet();
	}
	
	return output;
}

// Read NewVIKAR detector file and load bars into an array
// Returns the number of detectors loaded from the file
// Assumes the following detector file format for each bar in file
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
unsigned short ReadDetFile(const char* fname, Planar *bar_array){
	std::ifstream detfile(fname);
	if(!detfile.good()){ return 0; }
		
	std::vector<NewVIKARDet> detectors;
	std::string bar_type;
	float values[9];
	
	while(true){
		for(unsigned short i = 0; i < 6; i++){ 
			detfile >> values[i];
		}
		
		// Set the size of the bar
		detfile >> bar_type;
		if(!(bar_type == "small" || bar_type == "medium" || bar_type == "large")){
			for(unsigned short i = 6; i < 9; i++){ detfile >> values[i]; }
		}
		
		detectors.push_back(NewVIKARDet(values, bar_type));
		if(detfile.eof()){ break; }
	}	
	detfile.close();
	
	// Generate the Planar bar array
	bar_array = new Planar[detectors.size()];
	for(unsigned short i = 0; i < detectors.size(); i++){
		// Set the size of the bar
		if(detectors[i].bar_type == "small"){ bar_array[i].SetSmall(); }
		else if(detectors[i].bar_type == "medium"){ bar_array[i].SetMedium(); }
		else if(detectors[i].bar_type == "large"){ bar_array[i].SetLarge(); }
		else{ bar_array[i].SetSize(detectors[i].data[6],detectors[i].data[7],detectors[i].data[8]); }
		
		bar_array[i].SetPosition(detectors[i].data[0],detectors[i].data[1],detectors[i].data[2]); // Set the x,y,z position of the bar
		bar_array[i].SetBarRotation(detectors[i].data[3],detectors[i].data[4],detectors[i].data[5]); // Set the 3d rotation of the bar
	}
	
	return detectors.size();
}

// Perform a monte carlo simulation on an arbitrary configuration
// of VANDLE bars from an array. Returns the number of hits detected
// Generates two files...
//  xyz.dat - Contains 3-tuples of (x,y,z) for all detected hits
//  faces.dat - Contains 4-tuples of (bar#,face_x,face_y,face_z) for all detected hits
unsigned int TestDetSetup(Planar *bar_array, unsigned short num_bars, unsigned int num_trials){
	std::ofstream xyz("xyz.dat");
	std::ofstream faces("faces.dat");
	
	double tempx, tempy, tempz;
	unsigned int count;
	unsigned short bar;
	Vector3 temp_vector;
	Vector3 temp_ray;
	
	count = 0;
	for(unsigned int i = 0; i < num_trials; i++){
		UnitSphereRandom(temp_ray); // Generate a uniformly distributed random point on the unit sphere
		for(bar = 0; bar < num_bars; bar++){
			if(bar_array[bar].FaceIntersect(zero_vector, temp_ray, temp_vector, tempx, tempy, tempz) != -1){
				// A hit was detected on one of the faces
				xyz << temp_vector.axis[0] << "\t" << temp_vector.axis[1] << "\t" << temp_vector.axis[2] << "\n"; // Position data
				faces << bar << "\t" << tempx << "\t" << tempy << "\t" << tempz << "\n"; // Local face position data
				count++;
				break;
			}
		}
	}
	xyz.close();
	faces.close();
	return count;
}
