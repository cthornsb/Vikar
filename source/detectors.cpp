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
	material_id = 0;
	detX = Vector3(1.0, 0.0, 0.0);
	detY = Vector3(0.0, 1.0, 0.0);
	detZ = Vector3(0.0, 0.0, 1.0);
	length = 1.0;
	width = 1.0;
	depth = 1.0;
	small = false;
	med = false;
	large = false;
	need_set = true;
	use_recoil = false;
	use_eject = false;
	use_material = false;
	type = "unknown";
	subtype = "unknown";
}

// Set the global face coordinates (wrt global origin)
// Each vertex is the center coordinate of one of the faces 
void Planar::_set_face_coords(){
	// Calculate the center points of each face
	GlobalFace[0] = position-detZ*(depth/2.0); // Front face
	GlobalFace[1] = position+detX*(width/2.0); // Right face
	GlobalFace[2] = position+detZ*(depth/2.0); // Back face
	GlobalFace[3] = position-detX*(width/2.0); // Left face
	GlobalFace[4] = position+detY*(length/2.0); // Top face
	GlobalFace[5] = position-detY*(length/2.0); // Bottom face

	need_set = false;
}

// Return the unit vector of one of the faces. Returns the zero vector if the face is undefined
void Planar::GetUnitVector(unsigned int face_, Vector3 &unit){
	if(face_ == 0){ unit = detZ; } // +z local axis face (front)
	else if(face_ == 1){ unit = detX; } // +x local axis (right)
	else if(face_ == 2){ unit = detZ; unit *= -1; } // -z local axis (back)
	else if(face_ == 3){ unit = detX; unit *= -1; } // -x local axis (left)
	else if(face_ == 4){ unit = detY; } // +y local axis (top)
	else if(face_ == 5){ unit = detY; unit *= -1; } // -y local axis (bottom)
	else{ unit = Vector3(0.0, 0.0, 0.0); }
}	

// Return the local bar frame coordinates of a global coordinate
void Planar::GetLocalCoords(const Vector3 &world_coords, double &x, double &y, double &z){
	Vector3 temp = (world_coords - position);
	x = temp.Dot(detX);
	y = temp.Dot(detY);
	z = temp.Dot(detZ);
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
void Planar::SetRotation(double theta_, double phi_, double psi_){
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
}

// Manually set the local bar unit vectors
// I would advise against this, use SetRotation instead
// Note: This method does not update the bar rotation values
//  and should therefore only be used for testing and debugging
void Planar::SetUnitVectors(const Vector3 &unitX, const Vector3 &unitY, const Vector3 &unitZ){
	detX = unitX; detY = unitY; detZ = unitZ;
	detX.Normalize(); detY.Normalize(); detZ.Normalize();
	need_set = true;
}

// Check if a point (in local coordinates) is within the bounds of the primitive
// Return true if the coordinates are within the primitive and false otherwise
bool Planar::CheckBounds(unsigned int face_, double x_, double y_, double z_){
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

// Find if a ray (from the origin) intersects the infinite
// plane defined by one of the faces of the VANDLE bar
// face 0 is along the +z local axis
// face 1 is along the +x local axis
// face 2 is along the -z local axis
// face 3 is along the -x local axis
// face 4 is along the +y local axis
// face 5 is along the -y local axis
bool Planar::PlaneIntersect(const Vector3 &offset_, const Vector3 &direction_, unsigned int face_, Vector3 &P){
	if(need_set){ _set_face_coords(); }
	
	Vector3 unit;
	GetUnitVector(face_, unit);
	
	// The ray vector has the parametric form ray = offset_ + t*direction_
	// First find the intersection points between the ray and a plane containing the face polygon
	//double denom = direction_.Dot(unit);
	//if(denom < 1E-8){ return false; }
	
	double t = (GlobalFace[face_].Dot(unit)-offset_.Dot(unit))/(direction_.Dot(unit));
	P = offset_ + direction_*t; // The intersection point on the plane
	
	if(t >= 0){ return true; } // The ray intersects the plane
	return false; // The plane is behind the ray, the ray will never intersect it
}

// Determine if a ray from the origin intersected a face of the bar
// and return the face id number. Return -1 if no intersection was found
int Planar::FaceIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz){
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
int Planar::TopBotIntersect(const Vector3 &offset_, const Vector3 &direction_, Vector3 &intersect, double &px, double &py, double &pz){
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
bool Planar::IntersectPrimitive(const Vector3 &offset_, const Vector3 &direction_, Vector3 &P1, Vector3 &P2, 
								int &face1, int &face2, double &px, double &py, double &pz){
	if(need_set){ _set_face_coords(); }
	
	double px2, py2, pz2;
	int face_count = 0;
	Vector3 ray, unit;
	for(unsigned int i = 0; i < 6; i++){
		GetUnitVector(i, unit);
		if(PlaneIntersect(offset_, direction_, i, ray)){ // The ray intersects the plane containing this face
			// Transform the intersection point into local coordinates and check if they're within the bounds
			GetLocalCoords(ray, px2, py2, pz2);
			if(CheckBounds(i, px2, py2, pz2)){ // The face was struck
				if(face_count == 0){ 
					px = px2; py = py2; pz = pz2;
					P1 = ray; 
					face1 = i; 
				}
				else if(face_count == 1){ 
					P2 = ray; 
					face2 = i; 
					
					if(P2.Length() < P1.Length()){ px = px2; py = py2; pz = pz2; } // Ensure we get the hit on the plane facing the target
					//break;
				}
				face_count++;
			}
		} // if(PlaneIntersect(offset_, direction_, i, ray))
	} // for(unsigned int i = 0; i < 6; i++)
	
	if(face_count == 2){ return true; }
	else if(face_count == 0){ return false; }
	else{ std::cout << " Unexpected number of face intersects: " << face_count << std::endl; }
	return false; 
}

// Trace a ray through the detector and calculate the thickness it sees between two faces (f1 and f2)
// Return -1 if the ray does not travel through both faces
double Planar::GetApparentThickness(const Vector3 &offset_, const Vector3 &direction_, unsigned int f1_, unsigned int f2_, Vector3 &intersect1, Vector3 &intersect2){
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
	std::stringstream stream;
	for(unsigned int i = 0; i < 6; i++){
		stream << GlobalFace[i].axis[0] << "\t" << GlobalFace[i].axis[1] << "\t" << GlobalFace[i].axis[2] << "\n";
	}
	stream << position.axis[0] << "\t" << position.axis[1] << "\t" << position.axis[2];
	return stream.str();
}

// Dump VIKAR detector format string
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
std::string Planar::DumpDet(){
	if(need_set){ _set_face_coords(); }
	std::stringstream stream;
	stream << position.axis[0] << "\t" << position.axis[1] << "\t" << position.axis[2];
	stream << "\t" << theta << "\t" << phi << "\t" << psi;
	stream << "\t" << type << "\t" << subtype;
	stream << "\t" << length << "\t" << width << "\t" << depth;
	return stream.str();
}

/////////////////////////////////////////////////////////////////////
// Wall class
//  Used to generate flat walls of bars
/////////////////////////////////////////////////////////////////////

// Initialize the wall object.
// Spacing is the physical distance between the edges of bars.
// This method leaves a distance of spacing/2 on each side of the wall
// so that walls having the same bar spacing will align properly
void Wall::Initialize(unsigned int num_bars_, double spacing, double bar_length, double bar_width, double bar_depth){
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
	for(unsigned int i = 0; i < num_bars; i++){
		bars[i].SetSize(bar_length, bar_width, bar_depth);
		bars[i].SetRotation(theta, phi, psi);
		bars[i].SetPosition(position + detX*((-width/2.0+spacing/2.0+bar_width/2.0)+i*(bar_width+spacing)));
	}

	init = true;
}

// Return a pointer to a bar in this wall
// Return null if the bar does not exist
Planar* Wall::GetBar(unsigned int bar){ 
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
	for(unsigned int i = 0; i < num_bars; i++){
		output += bars[i].DumpDet();
	}
	
	return output;
}

// Read NewVIKAR detector file and load bars into an array
// Returns the number of detectors loaded from the file
// Assumes the following detector file format for each bar in file
// X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Bar_Type [Length(m) Width(m) Depth(m)]
unsigned int ReadDetFile(const char* fname_, std::vector<Planar*> &bar_vector){
	std::ifstream detfile(fname_);
	if(!detfile.good()){ return 0; }
	bar_vector.clear();
	
	std::vector<NewVIKARDet> detectors;
	std::string line;

	while(true){
		getline(detfile, line);
		if(detfile.eof()){ break; }
		if(line[0] == '#'){ continue; } // Commented line
	
		detectors.push_back(NewVIKARDet(line));
	}	
	detfile.close();
	
	Planar *current_det;
	for(std::vector<NewVIKARDet>::iterator iter = detectors.begin(); iter != detectors.end(); iter++){
		current_det = new Planar();
		if(iter->type == "vandle"){ // Vandle bars only detect ejectiles (neutrons)
			current_det->SetEjectile();
			if(iter->subtype == "small"){ current_det->SetSmall(); }
			else if(iter->subtype == "medium"){ current_det->SetMedium(); }
			else if(iter->subtype == "large"){ current_det->SetLarge(); }
			else{ current_det->SetSize(iter->data[6],iter->data[7],iter->data[8]); }
		}
		else{
			if(iter->type == "recoil"){ current_det->SetRecoil(); } // Recoil detectors only detect recoils
			else if(iter->type == "dual"){ // Dual detectors detect both recoils and ejectiles
				current_det->SetRecoil(); 
				current_det->SetEjectile();
			}
			current_det->SetSize(iter->data[6],iter->data[7],iter->data[8]);
		}
		
		// Set the position and rotation
		current_det->SetPosition(iter->data[0],iter->data[1],iter->data[2]); // Set the x,y,z position of the bar
		current_det->SetRotation(iter->data[3],iter->data[4],iter->data[5]); // Set the 3d rotation of the bar
		current_det->SetType(iter->type);
		current_det->SetSubtype(iter->subtype);
		bar_vector.push_back(current_det);
	}
	
	return bar_vector.size();
}

// Perform a monte carlo simulation on an arbitrary configuration
// of detectors from an array. Returns the number of hits detected
// Generates two files...
//  xyz.dat - Contains 3-tuples of (x,y,z) for all detected hits
//  faces.dat - Contains 3-tuples of (face_x,face_y,face_z) for all detected hits
unsigned int TestDetSetup(Planar *bar_array, unsigned int num_bars, unsigned int num_trials){
	std::ofstream xyz("xyz.dat");
	std::ofstream faces("faces.dat");
	
	double tempx, tempy, tempz;
	unsigned int count, total;
	unsigned int bar;
	Vector3 flight_path;
	Vector3 temp_vector1;
	Vector3 temp_vector2;
	Vector3 temp_ray;
	int face1, face2;
	
	double penetration, dist_traveled;
	double fpath1, fpath2;
		
	total = 0; count = 0;
	while(count < num_trials){
		UnitSphereRandom(temp_ray); // Generate a uniformly distributed random point on the unit sphere
		for(bar = 0; bar < num_bars; bar++){
			if(bar_array[bar].IntersectPrimitive(zero_vector, temp_ray, temp_vector1, temp_vector2, face1, face2, tempx, tempy, tempz)){
				if(bar_array[bar].IsEjectileDet()){
					flight_path = (temp_vector2-temp_vector1); // The vector pointing from the first intersection point to the second
					penetration = frand(); // The fraction of the bar which the neutron travels through
					dist_traveled = flight_path.Length()*penetration; // Random distance traveled through bar
					fpath1 = temp_vector1.Length(); // Distance from reaction to first intersection point
					fpath2 = temp_vector2.Length(); // Distance from reaction to second intersection point
			
					// Calculate the total distance traveled and the interaction point inside the detector
					if(fpath1 <= fpath2){ 
						dist_traveled += fpath1; 
						flight_path = temp_vector1 + flight_path*penetration;
					}
					else{ 
						dist_traveled += fpath2; 
						flight_path = temp_vector2 + flight_path*penetration;
					}
		
					xyz << flight_path.axis[0] << "\t" << flight_path.axis[1] << "\t" << flight_path.axis[2] << "\n"; // Position data
					faces << tempx << "\t" << tempy << "\t" << tempz << "\n"; // Local face position data
					count++;
				}
			}
		}
		total++;
	}
	xyz.close();
	faces.close();

	return total;
}
