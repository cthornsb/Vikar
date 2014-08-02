// Det_Planar.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:35:20 2014

#include "../include/planar.h"

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
	GlobalFace[0] = position-barZ*(depth/2.0); // Front face
	GlobalFace[1] = position+barX*(width/2.0); // Right face
	GlobalFace[2] = position+barZ*(depth/2.0); // Back face
	GlobalFace[3] = position-barX*(width/2.0); // Left face
	GlobalFace[4] = position+barY*(length/2.0); // Top face
	GlobalFace[5] = position-barY*(length/2.0); // Bottom face
	need_set = false;
}

// Return the local bar frame coordinates of a global coordinate
void Planar::GetLocalCoords(Vector3 world_coords, double &x, double &y, double &z){
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
void Planar::SetPosition(Vector3 pos){
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
void Planar::SetUnitVectors(Vector3 unitX, Vector3 unitY, Vector3 unitZ){
	barX = unitX; barY = unitY; barZ = unitZ;
	barX.Normalize(); barY.Normalize(); barZ.Normalize();
	need_set = true;
}

// Find if a ray (from the origin) intersects the infinite
// plane defined by one of the faces of the VANDLE bar
// face 0 is along the +z local axis
// face 1 is along the +x local axis
// face 2 is along the -z local axis
// face 3 is along the -x local axis
// face 4 is along the +y local axis
// face 5 is along the -y local axis
bool Planar::PlaneIntersect(Vector3 ray, unsigned short face, Vector3 &P_){
	if(need_set){ _set_face_coords(); }
	
	Vector3 unit;
	if(face == 0){ unit = barZ; } 			// +z local axis face (front)
	else if(face == 1){ unit = barX; } 		// +x local axis (right)
	else if(face == 2){ unit = barZ; unit *= -1; }  // -z local axis (back)
	else if(face == 3){ unit = barX; unit *= -1; }  // -x local axis (left)
	else if(face == 4){ unit = barY; } 		// +y local axis (top)
	else if(face == 5){ unit = barY; unit *= -1; }  // -y local axis (bottom)
	else{ return false; }
	
	double denom = unit.Dot(ray);
	if(std::abs(denom) > 1E-6){
		// Calculate the point of intersection
		double d = GlobalFace[face].Dot(unit) / denom;
		P_ = ray*d;
    		
		if(d >= 0){ return true; } 
		else{ return false; } // Negative d means the plane is behind us
	}
	return false;
}

// Determine if a ray from the origin intersected a face of the bar
// and return the face id number. Return -1 if no intersection was found
short Planar::FaceIntersect(Vector3 ray, Vector3 &intersect, double &px, double &py, double &pz){
	if(need_set){ _set_face_coords(); }
	
	// Check the Front face (Normal to local +z axis)
	if(PlaneIntersect(ray, 0, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (py >= -length/2.0 && py <= length/2.0)){ 
			return 0; // Ray intersects the front face				
		}
	}

	// Check the Right face (Normal to local +x axis)
	if(PlaneIntersect(ray, 1, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((pz >= -depth/2.0 && pz <= depth/2.0) && (py >= -length/2.0 && py <= length/2.0)){
			return 1; // Ray intersects the right face
		}
	}

	// Check the Back face (Normal to local -z axis)
	if(PlaneIntersect(ray, 2, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (py >= -length/2.0 && py <= length/2.0)){ 
			return 2; // Ray intersects the back face     				
		}
	}

	// Check the Left face (Normal to local -x axis)
	if(PlaneIntersect(ray, 3, intersect)){
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
short Planar::TopBotIntersect(Vector3 ray, Vector3 &intersect, double &px, double &py, double &pz){
	if(need_set){ _set_face_coords(); }
	
	// Check the Top face (Normal to local +y axis)
	if(PlaneIntersect(ray, 4, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (pz >= -depth/2.0 && pz <= depth/2.0)){ 
			return 4; // Ray intersects the top face				
		}
	}
	
	// Check the Bottom face (Normal to local -y axis)
	if(PlaneIntersect(ray, 5, intersect)){
		GetLocalCoords(intersect, px, py, pz);
		if((px >= -width/2.0 && px <= width/2.0) && (pz >= -depth/2.0 && pz <= depth/2.0)){ 
			return 5; // Ray intersects the bottom face				
		}
	}
	
	return -1;
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
		UnitRandom(temp_ray); // Generate a uniformly distributed random point on the unit sphere
		for(bar = 0; bar < num_bars; bar++){
			if(bar_array[bar].FaceIntersect(temp_ray, temp_vector, tempx, tempy, tempz) != -1){
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

/////////////////////////////////////////////////////////////////////
// DetAng_plan.f
/////////////////////////////////////////////////////////////////////

void DetAng_plan(double length, double width, unsigned short Nstrips, double z, double x, double y, 
		 bool xplane, bool yplane, double *xMax, double *xMin, double *yMax, double *yMin){ 
	// DetAng_plan 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for calculating the range of x and y spanned by each strip
	// of a strip detector planar coordinates (normal parallel to the z axis),
	// with the strip direction parallel to the x or y axis.
	// The length, width and strip division are passed from the parent program,
	// along with the z location of the plane, the xy coordinates of the detector
	// centre, and the direction of the strips.
	// The maximum and minumum x and y coordinates for each strip are calculated
	// and returned.
	
	// "---------------" << std::endl;
	// "Detector: " << detNo << " Strip: " << n << std::endl;
	// "DetXMax =" << "\t" << DetXMax[detNo][n] << std::endl;
	// "DetXMin =" << "\t" << DetXMin[detNo][n] << std::endl;
	// "DetYMax =" << "\t" << DetYMax[detNo][n] << std::endl;
	// "DetYMin =" << "\t" << DetYMin[detNo][n] << std::endl;
	
	unsigned short n;  
	for(n = 0; n < Nstrips; n++){
		// Calculate the position (cm) of the edges of the strip
		if(xplane){
			yMax[n] = (width/Nstrips)*n - width/2.0 + y; 
			yMin[n] = (width/Nstrips)*(n-1) - width/2.0 + y; 
			xMax[n] = x+(length/2.0); 
			xMin[n] = x-(length/2.0); 
		}
		else if(yplane){
			xMax[n] = (width/Nstrips)*n - width/2.0 + x; 
			xMin[n] = (width/Nstrips)*(n-1) - width/2.0 + x; 
			yMax[n] = y+(length/2.0); 
			yMin[n] = y-(length/2.0); 
		} 
	} 
} 

/////////////////////////////////////////////////////////////////////
// DetHit_plan.f
/////////////////////////////////////////////////////////////////////

void DetHit_plan(double **DetXMax, double **DetXMin, double **DetYMax, double **DetYMin, unsigned short Ndet_plan, unsigned short *Nstrips_plan,
		 double *DetZ_plan, double theta, double phi, unsigned short &hit, unsigned short &DetHit, unsigned short &StripHit){ 
	// DetHit_plan 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for determination of hits in planar detectors.
	// The precise range in x and y (cartesian coordinates)
	// spanned by the detector at the detector plane is calculated
	// to verify the hit.
	// The DetXMax,DetXMin,DetYMax,DetYMin are passed from the parent program,
	// along with the Ndet_plan, Nstrips_plan, z location and angle (theta,phi)
	// of the particle.
	// Hit, DetHit and StripHit are returned
	 
	double hit_r, hit_x, hit_y; 
	unsigned short n, s;  	 
	
	for(n = 0; n < Ndet_plan; n++){
		for(s = 0; s < Nstrips_plan[n]; s++){
			// For each detector, check the particle is heading towards it
			if (((theta*rad2deg) < 90.0 && DetZ_plan[n] > 0) || ((theta*rad2deg) > 90.0 && DetZ_plan[n] < 0)){
				hit_r = DetZ_plan[n]*std::tan(theta); // radius of hit from z-axis in detector plane
				hit_x = hit_r*std::sin(phi); // x position of hit in detector plane
				hit_y = hit_r*std::cos(phi); // y position of hit in detector plane
		
				// Check if particle is within detector area
				if (hit_x <= DetXMax[n][s] && hit_x >= DetXMin[n][s] && hit_y <= DetYMax[n][s] && hit_y >= DetYMin[n][s]){
					hit = hit + 1; 
					StripHit = s; 
					DetHit = n; 
				} 
			} // theta check
		} // Strip
	} // Detector
}

/////////////////////////////////////////////////////////////////////
// DetSet_plan_read.f
/////////////////////////////////////////////////////////////////////

void DetSet_plan_read(const char* fName, unsigned short Ndet, unsigned short *Nstrips, double *detLength, double *detWidth, double *DetZ, double *DetX, double *DetY,
		      bool *xPlane, bool *yPlane, double *dEthick, double *Ethick, double *Eres_dE, double *Pres_dE, double *Eres_E,
		      double *Pres_E, double *Pres_dE_energy, double *Pres_E_energy){ 
	// DetSet_plan_read 1.0 written by S.D.Pain on 7/6/2005
	//
	// Subroutine for reading in planar detector details and locations
	// from fName.
	// The planar detector coordinates are specified in cartesians,
	// and are read in in cm from the origin.
	// The input file format should be of the form:
	//
	// (Detector length (cm)) (Detector Width (cm))  (N. Strips)  (Z position (cm))
	// (X Position (cm)) (Y Position (cm))  (plane of strips ('x' or 'y')   (dE thickness (um))   (E thickness (um))
	// (dE Energy Resolution (MeV))  (dE Pos Resolution (cm))
	// (E Energy Resolution (MeV))  (E Pos Resolution (cm))
	// (dE Pos Resolution Energy (MeV))  (E Pos Resolution Energy (MeV))
	std::string plane;

	// DetSet_stat stores error status - T = good, F = bad
	std::ifstream file10(fName); 	
	
	// Read in the main data points from the SRIM output file
	for(unsigned short i = 0; i < Ndet; i++){
		file10 >> detLength[i] >> detWidth[i] >> Nstrips[i] >> DetZ[i];
		file10 >> DetX[i] >> DetY[i] >> plane >> dEthick[i] >> Ethick[i];
		file10 >> Eres_dE[i] >> Pres_dE[i] >> Eres_E[i] >> Pres_E[i];
		file10 >> Pres_dE_energy[i] >> Pres_E_energy[i]; 

		if(plane == "x"){
			xPlane[i] = true; 
			yPlane[i] = false;
		} 
		else if(plane == "y"){
			xPlane[i] = false; 
			yPlane[i] = true; 
		} 	
	} // Read in data loop
	
	file10.close(); 
} 

/////////////////////////////////////////////////////////////////////
// det_thick_plan.f
/////////////////////////////////////////////////////////////////////

void det_thick_plan(double theta, double phi, double DetThickness, double &thickness){ 
	// det_thick_plan 1.0 written by S.D.Pain on 8/6/2005
	//
	// Subroutine to calculate the detector thickness seen by a particle
	// for detectors in planar coordinates
	
	// Thickness due to theta
	thickness = DetThickness/(std::sqrt(pow(std::cos(theta), 2.0)) ); 
}

/////////////////////////////////////////////////////////////////////
// resolution_plan.f
/////////////////////////////////////////////////////////////////////

void resolution_plan(double energy, double theta, double phi, double Eres, double Pres, double PresE, double DetZ_plan, bool xPlane,
		     bool yPlane, double DetXMin, double DetYMin, double width, double length, unsigned short Nstrips, unsigned short StripHit, 
		     double newEnergy, double newTheta, double newPhi){ 
	// resolution_plan (VANDLE) 1.0 adapted by S.D.Pain on 2014/01/22
	// Subroutine to smear the particles by the detector resolution
	// in energy and position, for planar detectors
	
	// Pres = position resolution at 5.8MeV
	// Pres1 = position resolution at 'energy'
	
	double Pres1; 
	double dummy; 
	double hit_r, hit_x, hit_y; 
	
	double t_res; 
	double vel, n_mass, time, hit_vect; 
	
	// **** begin VANDLE specific ******
	n_mass = 1.0; 
	t_res = 3.0; 						// beam timing resolution in ns (Changed to 3 ns for vandle timing by Cory)
	vel = velocity(energy,n_mass); 
	// **** end VANDLE specific ******
	
	hit_r = DetZ_plan*tan(theta); 				// radius of hit from z-axis in detector plane
	hit_x = hit_r*sin(phi); 				// x position of hit in detector plane
	hit_y = hit_r*cos(phi); 				// y position of hit in detector plane
	
	// **** begin VANDLE specific ******
	hit_vect = std::sqrt(pow(DetZ_plan, 2)+pow(hit_r, 2)); 
	time = (((frand()-0.5)*5.0)+hit_vect)/vel; 		// 5.0 = bar thickness in cm (large bar)
	// time = (((rand(0)-0.5)*3.0)+hit_vect)/vel 		// 3.0 = bar thickness in cm (medium & small bar)
	time = time +rndgauss0(t_res); 
	// time = hit_r/vel 					// 2.5 = bar thickness in cm
	vel = hit_vect/time; 
	// **** end VANDLE specific ******
	
	// newEnergy = energy + rndgauss0(Eres) 		// removed for VANDLE
	newEnergy = 0.5*n_mass*pow(vel, 2 ); 			// added for VANDLE
	
	// Pres1 = (Pres/length*PresE)/energy*length 		// Pres at 'energy' !removed for VANDLE
	Pres1 = Pres;						// added for VANDLE
	
	if(xPlane){
		// hit_y = detyMin + (width/Nstrips*StripHit)-(width/Nstrips)/2.0	//assumes detyMin is bottom of detector, not strip
		hit_x = hit_x + rndgauss0(Pres1);
		hit_y = DetYMin + (width/Nstrips)/2.0; 
	}
	else if(yPlane){ 
		// hit_x = detxMin + (width/Nstrips*StripHit)-(width/Nstrips)/2.0 	// as above
		hit_x = DetXMin + (width/Nstrips)/2.0; 
		hit_y = hit_y + rndgauss0(Pres1); 
	} 
	
	Cart2Sphere(hit_x,hit_y,DetZ_plan,dummy,newTheta,newPhi); 
} 

/////////////////////////////////////////////////////////////////////
// strag_dE_plan.f
/////////////////////////////////////////////////////////////////////

void strag_dE_plan(double A, double Z, double detThick, double DetZ, double theta_old, double phi_old, double energy, 
		   double &theta_out, double &phi_out, double X, double conv_det){
	// strag_dE_cyl 1.0 written by S.D.Pain on 23/02/2005
	//
	// strag_dE_cyl 1.1 modified by S.D.Pain on 7/03/2005
	// to untilise transform2.0
	//
	// Subroutine for the calcualtion of the angular straggling in
	// the dE layer of a telescope, and application of its effect
	// on the measured angles in the E detector.
	//
	//     A = A of ion
	//     Z = Z of ion
	//     detThick = Apparent detector thickness in um (because of incident angle)
	//     DetPhi = azimuthal angle of detector
	//     theta_old = polar angle of incident ion (wrt target)
	//     phi_old = azimuthal angle of incident ion (wrt target)
	//     energy = energy of ion
	//     theta_new = polar angle of outgoing ion (wrt target)
	//     phi_new = azimuthal angle of outgoing ion (wrt target)
	//     X = radiation length of detector material
	//     conv_det = detector thickness conversion from um to mg/cm2
	//
	//     PhiDet = azimuthal angle of incidnent ion wrt detector

	double detThick_mg;
	double length1,length2;
	double theta_scat,phi_scat;
	double theta_scatW;
	double theta_new,phi_new;
	double oldX,oldY,oldZ;
	double newX,newY,newZ;
	double X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3;
	double DetSep,detAngle;
	double dummy,dummylength;

	// Calculate the detector thickness in mg/cm^2, and determine the
	// scattering width
	// detThick = 140 
	// energy = 5.0
	detThick_mg = detThick*conv_det;
	straggleA(theta_scatW,energy,Z,A,detThick_mg,X);
	theta_scat = rndgauss0(theta_scatW);
	theta_scat = sqrt(theta_scat*theta_scat*2.0);
	phi_scat = frand()*2.0*pi;
	dummylength = 1.0;

	// Direction of the particle
	// Convert the incident vector to cartesian, as this will be needed later
	Sphere2Cart(dummylength,theta_old,phi_old,oldX,oldY,oldZ);    

	// Transform the scattered ion vector to the lab coordinates (direction relative to target)    
	transform(theta_old,phi_old,theta_scat,phi_scat,theta_new,phi_new);
	     
	// Convert the scattered vector to cartesians. This will be used to
	// calculate the angle of the scattered ion relative to the detector. (direction again)
	Sphere2Cart(dummylength,theta_new,phi_new,newX,newY,newZ);

	// Separation of detectors (cm) possibly read this in at some pounsigned short 
	DetSep = 1.0;

	// Calculate absolute angle of the scattered ion wrt to detector plane 
	if(theta_new <= pi/2.0){ detAngle = theta_new; }
	else if(theta_new > pi/2.0){ detAngle = pi-theta_new; }

	// Calculate the length of the vector from the scattering pounsigned short to the
	// pounsigned short of incidence on the E detector
	length2 = DetSep/cos(detAngle);

	// Calculate the length of the vector from the target to
	// the scattering point
	length1 = DetZ/cos(theta_old);

	// the vector from target to scattering point
	unitV(oldX,oldY,oldZ,X1,Y1,Z1,dummy);
	X1 = X1*length1;
	Y1 = Y1*length1;
	Z1 = Z1*length1;

	// the vector from scattering pounsigned short to plane of E detector
	unitV(newX,newY,newZ,X2,Y2,Z2,dummy);
	X2 = X2*length2;
	Y2 = Y2*length2;
	Z2 = Z2*length2;

	// The vector from target to pounsigned short of incidence on E detector
	X3 = X1+X2;
	Y3 = Y1+Y2;
	Z3 = Z1+Z2;
	 
	// Convert this back to spherical polars, and return theta and phi
	Cart2Sphere(X3,Y3,Z3,dummylength,theta_out,phi_out);
}
