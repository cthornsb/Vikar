// detectors.cpp
// Cory Thornsberry

#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// NewVIKARdet
/////////////////////////////////////////////////////////////////////

NewVIKARdet::NewVIKARdet(){ 
	for(unsigned int i = 0; i < 9; i++){ data[i] = 0.0; }
	type = "unknown";
	subtype = "unknown"; 
	material = "none";
}

NewVIKARdet::NewVIKARdet(std::string input_){
	for(unsigned int i = 0; i < 9; i++){ data[i] = 0.0; }
	SetValues(input_);
}

void NewVIKARdet::SetValues(std::string input){
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

std::string NewVIKARdet::DumpDet(){
	std::stringstream stream;
	stream << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5];
	stream << "\t" << type << "\t" << subtype << "\t" << data[6] << "\t" << data[7] << "\t" << data[8];
	return stream.str();
}

/////////////////////////////////////////////////////////////////////
// Wall class
//  Used to generate flat walls of bars
/////////////////////////////////////////////////////////////////////

/** Initialize the wall object.
  * Spacing is the physical distance between the edges of bars.
  * This method leaves a distance of spacing/2 on each side of the wall
  * so that walls having the same bar spacing will align properly.
  */
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
	bars = new Primitive[num_bars];

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

/** Return a pointer to a bar in this wall
  * Return null if the bar does not exist.
  */
Primitive* Wall::GetBar(unsigned int bar){ 
	if(init && bar < num_bars){ return &bars[bar]; }
	return NULL;
}

/// Test the Primitive wall class.
void Wall::WallTest(){
	if(!init){ std::cout << " Initialized: No\n"; return; }
	else{ std::cout << " Initialized: Yes\n"; }
	
	std::cout << " Length: " << length << std::endl;
	std::cout << " Width: " << width << std::endl;
	std::cout << " Depth: " << depth << std::endl;
}

/** Dump VIKAR detector format string
  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
  */
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

/////////////////////////////////////////////////////////////////////
// Elliptical
/////////////////////////////////////////////////////////////////////

void Elliptical::SetRadii(const double &rLong_, const double &rShort_){
	rLong = (rLong_>=0.0)?rLong_:-1*rLong_;
	rShort = (rShort_>=0.0)?rShort_:-1*rShort_;
	length = rLong*2.0;
	width = rShort*2.0;
}

/** Check if a point (in local coordinates) is within the bounds of the primitive
  * Return true if the coordinates are within the primitive and false otherwise.
  * The Elliptical class checks if the face intersect is within 
  */
bool Elliptical::CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_){
	if(face_ == front_face || face_ == back_face){ // Front face or back face.
		if(x_*x_/(rLong*rLong) + y_*y_/(rShort*rShort) <= 1){ return true; }
	}
	return false;
}

/////////////////////////////////////////////////////////////////////
// Polygon
/////////////////////////////////////////////////////////////////////

/** Check if a point (in local coordinates) is within the bounds of the primitive
  * Return true if the coordinates are within the primitive and false otherwise.
  * The Elliptical class checks if the face intersect is within 
  */
bool Polygonal::CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_){
	if(face_ == front_face || face_ == back_face){ // Front face or back face.
		// Use the even-odd crossing test.
		if(poly.IsInside(x_, y_)){ return true; }
	}
	return false;
}

/////////////////////////////////////////////////////////////////////
// Annular
/////////////////////////////////////////////////////////////////////

void Annular::SetRadii(const double &inRadius_, const double &outRadius_){
	inRadius = (inRadius_>=0.0)?inRadius_:-1*inRadius_;
	outRadius = (outRadius_>=0.0)?outRadius_:-1*outRadius_;
	length = inRadius*2.0;
	width = outRadius*2.0;
}

/** Check if a point (in local coordinates) is within the bounds of the primitive
  * Return true if the coordinates are within the primitive and false otherwise.
  * The Elliptical class checks if the face intersect is within 
  */
bool Annular::CheckBounds(const unsigned int &face_, const double &x_, const double &y_, const double &z_){
	if(face_ == front_face || face_ == back_face){ // Front face or back face.
		double radius = std::sqrt(x_*x_ + y_*y_);
		//std::cout << inRadius << ", " << outRadius << ", " << radius << std::endl;
		if(radius >= inRadius && radius <= outRadius){ return true; }
	}
	return false;
}

/** Read NewVIKAR detector file and load detectors into a vector of pointers.
  * Returns the number of detectors loaded from the file.
  * Assumes the following detector file format for each detector in file
  * X(m) Y(m) Z(m) Theta(rad) Phi(rad) Psi(rad) Type Subtype Length(m) Width(m) Depth(m) Material.
  */
int ReadDetFile(const char* fname_, std::vector<Primitive*> &detectors){
	// Read VIKAR detector setup file or manually setup simple systems
	std::ifstream detfile(fname_);
	if(!detfile.good()){ return -1; }

	std::vector<NewVIKARdet*> temp_detectors;
	std::string line;

	while(true){
		getline(detfile, line);
		if(detfile.eof()){ break; }
		if(line[0] == '#'){ continue; } // Commented line

		temp_detectors.push_back(new NewVIKARdet(line));
	}	
	detfile.close();

	// Fill the detector
	for(std::vector<NewVIKARdet*>::iterator iter = temp_detectors.begin(); iter != temp_detectors.end(); iter++){
		if((*iter)->type == "vandle" || (*iter)->subtype == "planar"){ detectors.push_back(new Planar((*iter))); }
		else if((*iter)->subtype == "cylinder"){ detectors.push_back(new Cylindrical((*iter))); }
		else if((*iter)->subtype == "sphere"){ detectors.push_back(new Spherical((*iter))); }
		else if((*iter)->subtype == "ellipse"){ detectors.push_back(new Elliptical((*iter))); }
		else if((*iter)->subtype == "polygon"){ detectors.push_back(new Polygonal((*iter))); }
		else if((*iter)->subtype == "annular"){ detectors.push_back(new Annular((*iter))); }
		else{ 
			std::cout << " Unknown detector of type = " << (*iter)->type << " and subtype = " << (*iter)->subtype << std::endl;
			continue; 
		}
		delete (*iter);
	}
	
	return detectors.size();
}
