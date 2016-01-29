// vikar_lib.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:33:21 2014

#include "vikar_core.h"
#include "detectors.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals (for fortran commons)
/////////////////////////////////////////////////////////////////////

const double c = 2.99792458E8; // m/s
const double pi = 3.1415926540;
const double deg2rad = pi/180.0;
const double rad2deg = 180.0/pi;
const double LN2 = 0.6931471805;

/** This function finds the point in 2d space where two rays intersect.
  * Return the parameters t1 and t2 for the following parametric vector equations
  *  P1(t1) = p1_ + d1_*t1
  *  P2(t2) = p2_ + d2_*t2
  *  such that P1(t1) = P2(t2) = P
  */
bool GetT1T2(const Vector3 &p1_, const Vector3 &d1_, const Vector3 &p2_, const Vector3 &d2_, double &t1, double &t2){
	Vector3 V1 = p2_ - p1_;
	if(d1_.axis[0] != 0 && d1_.axis[1] != 0 && d2_.axis[0] != 0 && d2_.axis[1] != 0){ // General case.
		t2 = ((d1_.axis[0]/d1_.axis[1])*V1.axis[1] - V1.axis[0])*(1.0/(d2_.axis[0]*(1.0-d1_.axis[0]*d2_.axis[1]/(d1_.axis[1]*d2_.axis[0]))));
		t1 = (V1.axis[0]+d2_.axis[0]*t2)/d1_.axis[0];	
	}
	else if(d1_.axis[0] == 0){ // Line 1 is vertical.
		t2 = (p1_.axis[0] - p2_.axis[0]) / d2_.axis[0];
		t1 = (p2_.axis[1] - p1_.axis[1] + d2_.axis[1]*t2) / d1_.axis[1];
	}
	else if(d1_.axis[1] == 0){ // Line 1 is horizontal.
		t2 = (p1_.axis[1] - p2_.axis[1]) / d2_.axis[1];
		t1 = (p2_.axis[1] - p1_.axis[1] + d2_.axis[0]*t2) / d1_.axis[0];		
	}
	else if(d2_.axis[0] == 0){ // Line 2 is vertical.
		t1 = (p2_.axis[0] - p2_.axis[0]) / d1_.axis[0];
		t2 = (p1_.axis[1] - p2_.axis[1] + d1_.axis[1]*t1) / d2_.axis[1];
	}
	else if(d2_.axis[1] == 0){ // Line 2 is horizontal.
		t1 = (p2_.axis[1] - p2_.axis[1]) / d1_.axis[1];
		t2 = (p1_.axis[0] - p2_.axis[0] + d1_.axis[0]*t1) / d2_.axis[0];
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////
// Ray
/////////////////////////////////////////////////////////////////////

/** Construct a ray by supplying its starting point (x1_, y1_) and
  * a point through which it passes (x2_, y2_).
  */
Ray::Ray(const double &x1_, const double &y1_, const double &x2_, const double &y2_){
	pos = Vector3(x1_, y1_);
	Vector3 destination(x2_, y2_);
	dir = destination - pos;
}

/** Construct a ray by supplying its starting point pos_ and
  * and its direction (dx_, dy_).
  */
Ray::Ray(const Vector3 &pos_, const double &dx_, const double &dy_){
	pos = pos_;
	dir = Vector3(dx_, dy_);
}

/** Construct a ray by supplying its starting point (x_, y_) and
  * and its direction dir_.
  */
Ray::Ray(const double &x_, const double &y_, const Vector3 &dir_){
	pos = Vector3(x_, y_);
	dir = dir_;
}

/** Construct a ray by supplying its starting point pos_ and
  * and its direction dir_.
  */
Ray::Ray(const Vector3 &pos_, const Vector3 &dir_){
	pos = pos_;
	dir = dir_;
}

/// Construct a ray from a line segment.
Ray::Ray(const Line &line_){
	pos = line_.p1;
	dir = line_.p2 - line_.p1;
}

/// Assignment operator.
const Ray& Ray::operator = (const Ray &other){
	pos = other.pos;
	dir = other.dir;
	return *this;
}

/// Return true if this ray intersects another ray in 2d space.
bool Ray::Intersect(const Ray &other_, Vector3 &p){
	double t1, t2;
	if(GetT1T2(pos, dir, other_.pos, other_.dir, t1, t2)){
		p = pos + dir*t1;
		return (t1 >= 0.0 && t2 >= 0.0);
	}
	return false;
}

/// Return true if this ray intersects a line segment in 2d space.
bool Ray::Intersect(const Line &line_, Vector3 &p){
	double t1, t2;
	if(GetT1T2(pos, dir, line_.p1, line_.dir, t1, t2)){
		p = pos + dir*t1;
		return (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0));
	}
	return false;
}

/////////////////////////////////////////////////////////////////////
// Line
/////////////////////////////////////////////////////////////////////

/** Construct a line segment by supplying its starting point (x1_, y1_) and
  * its ending point (x2_, y2_).
  */
Line::Line(const double &x1_, const double &y1_, const double &x2_, const double &y2_){
	p1 = Vector3(x1_, y1_);
	p2 = Vector3(x2_, y2_);
	length = (p2 - p1).Length();
	dir = p2 - p1;
}

/** Construct a line segment by supplying its starting point pos_ its
  * direction (dx_, dy_) and its length_.
  */
Line::Line(const Vector3 &pos_, const double &dx_, const double &dy_, const double &length_){
	p1 = pos_;
	length = length_;
	p2 = p1 + Vector3(dx_, dy_)*length;
	dir = p2 - p1;
}

/** Construct a line segment by supplying its starting point (x_, y_) its
  * direction dir_ and its length_.
  */
Line::Line(const double &x_, const double &y_, const Vector3 &dir_, const double &length_){
	p1 = Vector3(x_, y_);
	length = length_;
	p2 = p1 + dir_ * length;
	dir = p2 - p1;
}

/** Construct a line segment by supplying its starting point pos_ its
  * direction dir_ and its length_.
  */
Line::Line(const Vector3 &pos_, const Vector3 &dir_, const double &length_){
	p1 = pos_;
	length = length_;
	p2 = p1 + dir_*length;
	dir = p2 - p1;
}

/// Construct a line segment from a ray by specifying its length_.
Line::Line(const Ray &ray_, const double &length_){
	p1 = ray_.pos;
	dir = ray_.dir;
	length = length_;
	p2 = p1 + dir*length;
}

/// Assignment operator.
const Line& Line::operator = (const Line& other_){
	p1 = other_.p1;
	p2 = other_.p2;
	dir = other_.dir;
	length = other_.length;
	return *this;
}

/// Return true if this line segment intersects a ray in 2d space.
bool Line::Intersect(const Line &other_, Vector3 &p){
	double t1, t2;
	if(GetT1T2(p1, dir, other_.p1, other_.dir, t1, t2)){
		p = p1 + dir*t1;
		return ((t1 >= 0.0 && t1 <= 1.0) && (t2 >= 0.0 && t2 <= 1.0));
	}
	return false;
}

/// Return true if this line segment intersects another line segment in 2d space.
bool Line::Intersect(const Ray &ray_, Vector3 &p){
	double t1, t2;
	if(GetT1T2(p1, dir, ray_.pos, ray_.dir, t1, t2)){
		p = p1 + dir*t1;
		return ((t1 >= 0.0 && t1 <= 1.0) && t2 >= 0.0);
	}
	return false;
}

/////////////////////////////////////////////////////////////////////
// Vector3
/////////////////////////////////////////////////////////////////////

const Vector3& Vector3::operator = (const Vector3 &other){
	axis[0] = other.axis[0];
	axis[1] = other.axis[1];
	axis[2] = other.axis[2];
	return *this;
}

// Vector addition
const Vector3& Vector3::operator += (const Vector3 &other){
	axis[0] += other.axis[0];
	axis[1] += other.axis[1];
	axis[2] += other.axis[2];
	return *this;
}

// Vector subtraction
const Vector3& Vector3::operator -= (const Vector3 &other){
	axis[0] -= other.axis[0];
	axis[1] -= other.axis[1];
	axis[2] -= other.axis[2];
	return *this;
}

// Scalar multiplication
const Vector3& Vector3::operator *= (const double &scalar){
	axis[0] *= scalar;
	axis[1] *= scalar;
	axis[2] *= scalar;
	return *this;
}

Vector3 Vector3::operator + (const Vector3 &other) const {
	return Vector3(axis[0]+other.axis[0], axis[1]+other.axis[1], axis[2]+other.axis[2]);
}

Vector3 Vector3::operator - (const Vector3 &other) const {
	return Vector3(axis[0]-other.axis[0], axis[1]-other.axis[1], axis[2]-other.axis[2]);
}

Vector3 Vector3::operator * (const double &scalar) const {
	return Vector3(axis[0]*scalar, axis[1]*scalar, axis[2]*scalar);
}

// Dot product
double Vector3::Dot(const Vector3 &other) const {
	return (axis[0]*other.axis[0] + axis[1]*other.axis[1] + axis[2]*other.axis[2]);
}

// Return the cosine of the angle between two vectors.
double Vector3::CosAngle(const Vector3 &other) const {
	return (this->Dot(other)/(this->Length()*other.Length()));
}

// Cross product
Vector3 Vector3::Cross(const Vector3 &other) const {
	return Vector3((axis[1]*other.axis[2]-other.axis[1]*axis[2]),
		       (other.axis[0]*axis[2]-axis[0]*other.axis[2]),
		       (axis[0]*other.axis[1]-other.axis[0]*axis[1]));
}

// Return the length of the vector
double Vector3::Length() const {
	return std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
}

// Return the square of the vector
double Vector3::Square() const {
	return axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
}

double Vector3::Distance(const Vector3 &other) const {
	double x = axis[0]-other.axis[0];
	double y = axis[1]-other.axis[1];
	double z = axis[2]-other.axis[2];
	return std::sqrt(x*x+y*y+z*z);
}

// Normalize the vector and return the normalization parameter
double Vector3::Normalize(){
	double parameter = Length();
	axis[0] = axis[0]/parameter;
	axis[1] = axis[1]/parameter;
	axis[2] = axis[2]/parameter;
	return parameter;
}

std::string Vector3::Dump() const { 
	std::stringstream stream;
	stream << axis[0] << ", " << axis[1] << ", " << axis[2]; 
	return stream.str();
}

/////////////////////////////////////////////////////////////////////
// Matrix3
/////////////////////////////////////////////////////////////////////

void Matrix3::_initialize(){
	for(unsigned int i = 0; i < 3; i++){ 
		components[i][0] = 0.0; 
		components[i][1] = 0.0; 
		components[i][2] = 0.0; 
	}
}

Matrix3::Matrix3(){ _initialize(); }

Matrix3::Matrix3(double theta_, double phi_){
	_initialize();
	SetRotationMatrixSphere(theta_, phi_);
}

Matrix3::Matrix3(const Vector3 &vector_){
	_initialize();
	SetRotationMatrixSphere(vector_);
}

void Matrix3::SetRotationMatrixSphere(double theta_, double phi_, double psi_/*=0.0*/){
	double sin_theta = std::sin(theta_), cos_theta = std::cos(theta_);
	double sin_phi = std::sin(phi_), cos_phi = std::cos(phi_);
	double sin_psi = std::sin(psi_), cos_psi = std::cos(psi_);
	
	// Pitch-Roll-Yaw convention
	// Rotate by angle theta about the y-axis
	//  angle phi about the z-axis
	//  angle psi about the x-axis
	SetRow1(cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta); // Width axis
	SetRow2(sin_psi*sin_theta*cos_phi-cos_psi*sin_phi, sin_psi*sin_theta*sin_phi+cos_psi*cos_phi, cos_theta*sin_psi); // Length axis
	SetRow3(cos_psi*sin_theta*cos_phi+sin_psi*sin_phi, cos_psi*sin_theta*sin_phi-sin_psi*cos_phi, cos_theta*cos_psi); // Depth axis
}

void Matrix3::SetRotationMatrixSphere(const Vector3 &vector_){
	SetRotationMatrixSphere(vector_.axis[1], vector_.axis[2]);
}

void Matrix3::SetRotationMatrixCart(double x_, double y_, double z_){
	double r, theta, phi;
	Cart2Sphere(x_, y_, z_, r, theta, phi);
	SetRotationMatrixSphere(theta, phi);
}

void Matrix3::SetRotationMatrixCart(const Vector3 &vector_){
	SetRotationMatrixCart(vector_.axis[0], vector_.axis[1], vector_.axis[2]);
}

// Transform an input vector by this matrix
// Note: Expects the input vector to be in cartesian coordinates
void Matrix3::Transform(Vector3 &vector){
	double x = vector.axis[0], y = vector.axis[1], z = vector.axis[2];
	vector.axis[0] = components[0][0]*x + components[0][1]*y + components[0][2]*z;
	vector.axis[1] = components[1][0]*x + components[1][1]*y + components[1][2]*z;
	vector.axis[2] = components[2][0]*x + components[2][1]*y + components[2][2]*z;
}

// Transform an input vector by the transpose of this matrix
// Note: Expects the input vector to be in cartesian coordinates
void Matrix3::Transpose(Vector3 &vector){
	double x = vector.axis[0], y = vector.axis[1], z = vector.axis[2];
	vector.axis[0] = components[0][0]*x + components[1][0]*y + components[2][0]*z;
	vector.axis[1] = components[0][1]*x + components[1][1]*y + components[2][1]*z;
	vector.axis[2] = components[0][2]*x + components[1][2]*y + components[2][2]*z;
}

void Matrix3::Dump(){
	std::cout << " [" << components[0][0] << "\t" << components[0][1] << "\t" << components[0][2] << "]\n";
	std::cout << " [" << components[1][0] << "\t" << components[1][1] << "\t" << components[1][2] << "]\n";
	std::cout << " [" << components[2][0] << "\t" << components[2][1] << "\t" << components[2][2] << "]\n";
}

/////////////////////////////////////////////////////////////////////
// AngularDist
/////////////////////////////////////////////////////////////////////

/// Default constructor.
AngularDist::AngularDist(){
	reaction_xsection = 0.0;
	num_points = 0;
	init = false;
	
	com_theta = NULL;
	dsigma_domega = NULL;
	integral = NULL;
}

/// Destructor.
AngularDist::~AngularDist(){
	if(com_theta){ delete[] com_theta; }
	if(dsigma_domega){ delete[] dsigma_domega; }
	if(integral){ delete[] integral; }
}

/** Setup the angular distribution by reading it from a file.
  * Return true if the file is correctly loaded and contains a non-zero number
  * of data points and return false otherwise.
  * param[in] fname Filename of the angular distribution file. File should contain
  *  two columns. First is the center of mass angle (in degrees) and second is the
  *  differential cross section (in mb/Sr) at that CoM angle.
  * param[in] mtarg The molar mass of the target (g/mol).
  * param[in] tgt_thickness The thickness of the target (mg/cm^2).
  * param[in] beam_current The intensity of the beam (pps).
  */
bool AngularDist::Initialize(const char* fname, const double &mtarg, const double &tgt_thickness, const double &beam_current){
	if(init){ return false; }
	
	std::ifstream inFile(fname);
	if(!inFile.good()){ return false; }
	
	double x, y;
	std::vector<double> xvec, yvec;
	while(true){
		inFile >> x >> y;
		if(inFile.eof()){ break; }
		xvec.push_back(x);
		yvec.push_back(y);
		num_points++;
	}
	
	inFile.close();
	
	// Need at least two points to calculate total reaction X-section
	if(num_points > 1){		
		com_theta = new double[num_points];
		dsigma_domega = new double[num_points];
		integral = new double[num_points];

		unsigned int index = 0;
		std::vector<double>::iterator iter1, iter2;
		for(iter1 = xvec.begin(), iter2 = yvec.begin(); iter1 != xvec.end() && iter2 != yvec.end(); iter1++, iter2++){
			com_theta[index] = *iter1*deg2rad;
			dsigma_domega[index] = *iter2;
			index++;
		}

		reaction_xsection = 0.0;
		integral[0] = 0.0;
		
		// Calculate the reaction cross-section from the differential cross section
		double x1, x2, y1, y2;
		for(unsigned int i = 0; i < num_points-1; i++){
			x1 = com_theta[i]*deg2rad; y1 = dsigma_domega[i]*std::sin(x1);
			x2 = com_theta[i+1]*deg2rad; y2 = dsigma_domega[i+1]*std::sin(x2);
			reaction_xsection += 0.5*(x2-x1)*(y2+y1);
			integral[i+1] = reaction_xsection*2*pi; // The cumulative integral
		}
		
		reaction_xsection *= 2*pi;
		//rate = avagadro*tgt_thickness*reaction_xsection*(1E-27)/(500*mtarg); // Reaction probability
		rate = 0.0;
		rate *= beam_current; // Reaction rate (pps)
		return (init = true);
	}
	
	return false;
}

/** Setup the angular distribution using arrays.
  * Return true if the file is correctly loaded and contains a non-zero number
  * of data points and return false otherwise.
  * param[in] num_points_ The number of points in the input arrays.
  * param[in] angle_ Pointer to the array of center of mass angles (rad).
  * param[in] xsection_ Pointer to the array of differential cross sections (mb/Sr).
  */
bool AngularDist::Initialize(const unsigned int &num_points_, double *angle_, double *xsection_){
	if(init || num_points_ <= 1){ return false; }
	
	num_points = num_points_;

	com_theta = new double[num_points];
	dsigma_domega = new double[num_points];
	integral = new double[num_points];
	
	reaction_xsection = 0.0;
	integral[0] = 0.0;
	
	// Load the cross section values into the arrays
	for(unsigned int i = 0; i < num_points; i++){
		com_theta[i] = angle_[i]*deg2rad;
		dsigma_domega[i] = xsection_[i];
	}
	
	// Calculate the reaction cross-section from the differential cross section
	double x1, x2, y1, y2;
	for(unsigned int i = 0; i < num_points-1; i++){
		x1 = com_theta[i]*deg2rad; y1 = dsigma_domega[i]*std::sin(x1);
		x2 = com_theta[i+1]*deg2rad; y2 = dsigma_domega[i+1]*std::sin(x2);
		reaction_xsection += 0.5*(x2-x1)*(y2+y1);
		integral[i+1] = reaction_xsection*2*pi; // The cumulative integral
	}
	
	reaction_xsection *= 2*pi;
	
	return (init = true);
}

/// Setup the angular distribution using an isotropic distribution.
bool AngularDist::Initialize(const double &xsection_){
	if(init || xsection_ <= 0.0){ return false; }
	
	reaction_xsection = xsection_;
	
	return (init = true);	
}

/** Return a random angle sampled from the distribution (rad).
  * return the center of mass angle (rad) sampled from the distributionj
  * and return -1 if the sampling fails for any reason.
  */
double AngularDist::Sample(){
	if(!init){ return -1; }
	
	if(num_points > 0){ // Standard (non-isotropic) cross section.
		double rand_xsect = frand()*reaction_xsection;
		for(unsigned int i = 0; i < num_points-1; i++){
			if(integral[i] <= rand_xsect && rand_xsect <= integral[i+1]){ 
				return (com_theta[i] + (rand_xsect-integral[i])*(com_theta[i+1]-com_theta[i])/(integral[i+1]-integral[i]));
			}
		}
	}
	else{ // Isotropic cross section.
		return (frand()*pi);
	}
	
	return -1;
}

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

// Return true if an input string is in a vector of strings and false otherwise.
bool IsInVector(const std::string &input_, const std::vector<std::string> &str_vector_){
	for(std::vector<std::string>::const_iterator iter = str_vector_.begin(); iter != str_vector_.end(); iter++){
		if(input_ == *iter){ return true; }
	}
	return false;
}

// Get a random point on a gaussian beam profile
// fwhm_ is the FWHM of the beamspot in m
// offset_ is the offset in the negative z-direction (in m)
// beam is a 2d vector in the xy-plane (z=0) pointing from the origin to a point inside the target beamspot
void RandomGauss(double fwhm_, double offset_, Vector3 &beam){
	// Uniformly sample the gaussian profile
	double ranX = rndgauss0(fwhm_);
	double ranY = rndgauss0(fwhm_);
	
	beam = Vector3(ranX, ranY, -offset_);
}

// Get a random point on a circular beam profile
// radius_ is the beamspot radius in m
// offset_ is the offset in the negative z-direction (in m)
// beam is a 2d vector in the xy-plane (z=0) pointing from the origin to a point inside the target beamspot
void RandomCircle(double radius_, double offset_, Vector3 &beam){
	// Uniformly sample the circular profile
	double ranT = 2*pi*frand();
	double ranU = frand() + frand();
	double ranR;
	
	if(ranU > 1){ ranR = 2 - ranU; }
	else{ ranR = ranU; }
	ranR *= radius_;
	
	beam = Vector3(ranR*std::cos(ranT), ranR*std::sin(ranT), -offset_);
}

// Get a random point along the perimeter of a circle
// radius_ is the beamspot radius in m
// offset_ is the offset in the negative z-direction (in m)
// beam is a 2d vector in the xy-plane (z=0) pointing from the origin to a point inside the target beamspot
void RandomHalo(double radius_, double offset_, Vector3 &beam){
	// Uniformly sample the circular profile
	double ranT = 2*pi*frand();
	
	beam = Vector3(radius_*std::cos(ranT), radius_*std::sin(ranT), -offset_);
}

bool SetBool(std::string input_, std::string text_, bool &output){
	int idummy = atoi(input_.c_str());
	if(idummy == 1){ output = true; }
	else{ output = false; }
	std::cout << text_;
	if(output){ std::cout << ": Yes\n"; }
	else{ std::cout << ": No\n"; }
	return output;
}

bool SetBool(std::string input_, bool &output){
	int idummy = atoi(input_.c_str());
	if(idummy == 1){ output = true; }
	else{ output = false; }
	return output;
}

bool Prompt(std::string prompt_){
	std::string temp_input;
	while(true){
		std::cout << prompt_ << " (yes/no) "; std::cin >> temp_input;
		if(temp_input == "yes" || temp_input == "y"){ return true;; }
		else if(temp_input == "no" || temp_input == "n"){ return false; }
		else{ std::cout << "  Type yes or no\n"; }
	}
}

// Parse an input string and return text up to the first 
// occurance of white space or the first occurance of a '#'
std::string Parse(std::string input){
	std::string output = "";
	for(unsigned int i = 0; i < input.size(); i++){
		if(input[i] == ' ' || input[i] == '\t' || input[i] == '#'){ break; }
		output += input[i];
	}
	return output;
}

// Find the order of a number in powers of 10
double Order(double input_){
	double output = 1.0;
	if(input_ > 1.0){ // Of order 1.0E+x
		while(input_/output > 1.0){ output *= 10.0; }
	}
	else{ // Of order 1.0E-x
		while(input_/output < 1.0){ output = output/10.0; }
	}
	return output;
}

// Double absolute value
double dabs(double value){
	if(value < 0){ return -1.0*value; }
	return value;
}

// Return the minimum value
double min(double v1, double v2){
	if(v1 <= v2){ return v1; }
	else{ return v2; }
}

// Return the maximum value
double max(double v1, double v2){
	if(v1 > v2){ return v1; }
	else{ return v2; }
}

// Return a random number between low and high
double frand(double low, double high){
	return low+(double(rand())/RAND_MAX)*(high-low);
}

// Mimic the fortran rand() function
double frand(){
	return double(rand())/RAND_MAX;
}

// Sample a point on the surface of the unit sphere
void UnitSphereRandom(Vector3 &vec){
	double u = 2*frand()-1;
	double theta = 2*pi*frand();
	vec.axis[0] = std::sqrt(1-u*u)*std::cos(theta);
	vec.axis[1] = std::sqrt(1-u*u)*std::sin(theta);
	vec.axis[2] = u;
}

// Sample a point on the surface of the unit sphere
void UnitSphereRandom(double &theta, double &phi){
	phi = 2*pi*frand();
	theta = std::acos(2*frand()-1);
}

// Sample a point on the unit circle
double UnitCircleRandom(){
	return 2*pi*frand();
}

// Calculate proper bar spacing for a wall of VANDLE bars
// Leave half gaps at either edge for clearance to other walls
double BarSpacing(double total_width, double bar_width, unsigned int num_bars){
	return (total_width-num_bars*bar_width)/num_bars;
}

// Calculate the angular spacing between adjacent bars
double BarSpacingAngle(double radius, double spacing){
	return 2*std::asin(spacing/(2*radius));
}

// Wrap a value between min_val and max_val
double WrapValue(double value, double min_val, double max_val){
	if(value < min_val){ return max_val-(min_val-value); }
	else if(value > max_val){ return min_val+(value-max_val); }
	else{ return value; }
}

// Return the number of lines in a file
unsigned int GetLines(const char* input){
	std::ifstream file(input); 
	unsigned int count = 0;
	std::string line;
	while(std::getline(file,line)){ 
		count++; 
	}
	file.close();
	
	return count;
}

// Linearly interpolate between points
// Return the value y = f(x)
double Interpolate(double x1, double y1, double x2, double y2, double x){
	return ((y2-y1)/(x2-x1))*(x-x1)+y1;
}

// Return the distance between two points in 3d space
double Dist3d(const Vector3 &v1, const Vector3 &v2){
	return (v2-v1).Length();
}

/////////////////////////////////////////////////////////////////////
// Cart2Sphere.f
/////////////////////////////////////////////////////////////////////

void Cart2Sphere(double x, double y, double z, double &r, double &theta, double &phi){ 
	// Cart2Sphere 1.0 written by S.D.Pain on 4/12/2004
	//
	// Subroutine for converting a vector from cartesian coordinates
	// to spherical polar coordinates
	//
	// (x,y,z) are passed in, and (r,theata,phi) are calculated and
	// returned
	if(x == 0.0 && y == 0.0 && z == 0.0){
		r = 0.0; theta = 0.0; phi = 0.0;
	}
	else{
		r = std::sqrt(x*x + y*y + z*z);
		theta = std::acos(z/r);
		
		if(x == 0.0 && y == 0.0){ phi = 0.0; }
		else{ 
			double temp = std::sqrt(x*x + y*y); 
			if(x >= 0.0){ phi = std::acos(y/temp); }
			else{ phi = 2.0*pi - std::acos(y/temp); }
		}
	}
} 

void Cart2Sphere(double x, double y, double z, Vector3 &sphere){
	if(x == 0.0 && y == 0.0 && z == 0.0){
		sphere.axis[0] = 0.0; sphere.axis[1] = 0.0; sphere.axis[2] = 0.0;
	}
	else{
		sphere.axis[0] = std::sqrt(x*x + y*y + z*z);
		sphere.axis[1] = std::acos(z/sphere.axis[0]);
		
		if(x == 0.0 && y == 0.0){ sphere.axis[2] = 0.0; }
		else{ 
			double temp = std::sqrt(x*x + y*y); 
			if(x >= 0.0){ sphere.axis[2] = std::acos(y/temp); }
			else{ sphere.axis[2] = 2.0*pi - std::acos(y/temp); }
		}
	}
}

void Cart2Sphere(const Vector3 &cart, double &r, double &theta, double &phi){
	if(cart.axis[0] == 0.0 && cart.axis[1] == 0.0 && cart.axis[2] == 0.0){
		r = 0.0; theta = 0.0; phi = 0.0;
	}
	else{
		r = std::sqrt(cart.axis[0]*cart.axis[0] + cart.axis[1]*cart.axis[1] + cart.axis[2]*cart.axis[2]);
		theta = std::acos(cart.axis[2]/r);
		
		if(cart.axis[0] == 0.0 && cart.axis[1] == 0.0){ phi = 0.0; }
		else{ 
			double temp = std::sqrt(cart.axis[0]*cart.axis[0] + cart.axis[1]*cart.axis[1]); 
			if(cart.axis[0] >= 0.0){ phi = std::acos(cart.axis[1]/temp); }
			else{ phi = 2.0*pi - std::acos(cart.axis[1]/temp); }
		}
	}
}

void Cart2Sphere(const Vector3 &cart, Vector3 &sphere){
	if(cart.axis[0] == 0.0 && cart.axis[1] == 0.0 && cart.axis[2] == 0.0){
		sphere.axis[0] = 0.0; sphere.axis[1] = 0.0; sphere.axis[2] = 0.0;
	}
	else{
		sphere.axis[0] = std::sqrt(cart.axis[0]*cart.axis[0] + cart.axis[1]*cart.axis[1] + cart.axis[2]*cart.axis[2]);
		sphere.axis[1] = std::acos(cart.axis[2]/sphere.axis[0]);
		
		if(cart.axis[0] == 0.0 && cart.axis[1] == 0.0){ sphere.axis[2] = 0.0; }
		else{ 
			double temp = std::sqrt(cart.axis[0]*cart.axis[0] + cart.axis[1]*cart.axis[1]); 
			if(cart.axis[0] >= 0.0){ sphere.axis[2] = std::acos(cart.axis[1]/temp); }
			else{ sphere.axis[2] = 2.0*pi - std::acos(cart.axis[1]/temp); }
		}
	}
}

void Cart2Sphere(Vector3 &sphere){
	Vector3 cart = sphere;
	Cart2Sphere(cart, sphere);
}

/////////////////////////////////////////////////////////////////////
// Sphere2Cart.f
/////////////////////////////////////////////////////////////////////

void Sphere2Cart(double r, double theta, double phi, double &x, double &y, double &z){ 
	// Sphere2Cart 1.0 written by S.D.Pain on 4/12/2004
	//
	// Subroutine for converting a vector from spherical polar coordinates
	// to cartesian coordinates
	//
	// (r,theata,phi) are passed in, and (x,y,z) are calculated and
	// returned
	x = r*std::sin(theta)*std::cos(phi); 
	y = r*std::sin(theta)*std::sin(phi); 
	z = r*std::cos(theta);
}

void Sphere2Cart(double r, double theta, double phi, Vector3 &cart){
	cart.axis[0] = r*std::sin(theta)*std::cos(phi);
	cart.axis[1] = r*std::sin(theta)*std::sin(phi); 
	cart.axis[2] = r*std::cos(theta);
}

void Sphere2Cart(const Vector3 &sphere, double &x, double &y, double &z){
	x = sphere.axis[0]*std::sin(sphere.axis[1])*std::cos(sphere.axis[2]);
	y = sphere.axis[0]*std::sin(sphere.axis[1])*std::sin(sphere.axis[2]); 
	z = sphere.axis[0]*std::cos(sphere.axis[1]);
}

void Sphere2Cart(const Vector3 &sphere, Vector3 &cart){
	cart.axis[0] = sphere.axis[0]*std::sin(sphere.axis[1])*std::cos(sphere.axis[2]);
	cart.axis[1] = sphere.axis[0]*std::sin(sphere.axis[1])*std::sin(sphere.axis[2]); 
	cart.axis[2] = sphere.axis[0]*std::cos(sphere.axis[1]);
}

void Sphere2Cart(Vector3 &cart){
	Vector3 sphere = cart;
	Sphere2Cart(sphere, cart);
}

/////////////////////////////////////////////////////////////////////
// radlength.f
/////////////////////////////////////////////////////////////////////

double radlength(unsigned int A, unsigned int Z){
	// radlength 1.0 written by S.D.Pain on 11/02/2005
	// Function to calculate the radiation length of a material
	// in mg/cm^2
	// See Barnett et al., Phys. Rev. D 54 (1996) 1, page 135
	
	return 7.164e5*A/(Z*(Z+1.0)*std::log(287.0/(std::sqrt(Z)))); 
} 

/////////////////////////////////////////////////////////////////////
// rndgauss.f
/////////////////////////////////////////////////////////////////////

double rndgauss0(double w){
	// rndgauss0, a cut down version of rndgauss1; 
	// returns a random number with FWHM w centred at 0;

	if(w <= 0.0){ return 0; }
	
	double t = 0.0, tsq = 0.0;  
	const double c0=2.515517, c1=0.802853, c2=0.010328; 
	const double d1=1.432788, d2=0.189269, d3=0.001308; 
	const double widthfact=0.424628450; 

	t = frand(); 

	if(t > 0.5){ t = t-0.5; }
	if(t < 1e-30){ t = 11.46380587; }
	else{ 
		tsq = -log(t*t); 
		t = std::sqrt(tsq); 
	} 
	
	//     compute inverse by equn 26.2.23 
	t=t-(c0+c1*t+c2*tsq)/(1.00+d1*t+(d2+d3*t)*tsq);
	
	//     now randomize x in positive and negative direction
	if(frand() > 0.5){ t = -t; }
	
	//     now correct for standard deviation
	return widthfact*w*t; 
} 

/////////////////////////////////////////////////////////////////////
// straggleA.f
/////////////////////////////////////////////////////////////////////

void straggleA(double &theta, double energy, double Z, double A, double thickness, double X){ 
	// straggleA 1.0 written by S.D.Pain on 20/01/2004
	//
	// Subroutine to calculate the width of a gaussian distribution
	// of angles from the straggling of an energetic ion in a medium
	//   theta = sigma of distribution (spatial)
	//   Energy = energy of particle
	//   thickness = thickness of material
	//   X = radiation length of stopping material
	//   A = Mass number of ion
	//   Z = charge of ion
	
	// CURRENTLY ONLY TESTED FOR A LIMITED RANGE OF IONS, ENERGIES and TARGETS

	theta = 13.6/(std::sqrt(2.0*energy/A)*std::sqrt(2.0*energy*A))*Z*std::sqrt(thickness/X)*(1.0+0.038*std::log(thickness/X)); 
	theta = theta*std::sqrt(2.0); 
} 
