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

/////////////////////////////////////////////////////////////////////
// Vector3 Struct
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
// Matrix3 Struct
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

void Matrix3::SetRotationMatrixSphere(double theta_, double phi_){
	double sin_theta = std::sin(theta_), cos_theta = std::cos(theta_);
	double sin_phi = std::sin(phi_), cos_phi = std::cos(phi_);
	
	// Rz(phi)Ry(theta) rotation matrix
	SetRow1(cos_phi*cos_theta, -sin_phi, cos_phi*sin_theta);
	SetRow2(sin_phi*cos_theta, cos_phi, sin_phi*sin_theta);
	SetRow3(-sin_theta, 0.0, cos_theta);
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

// Load the differential cross section from a file.
// Return true if the file is correctly loaded and contains a non-zero number
// of data points. Does nothing if AngularDist has already been initialized
// Expects CoM angle in degrees and differential cross section in mb/sr
//
// Thickness (tgt_thickness) must be in units of mg/cm^2
// Current (beam_current) must be in units of particles/second (pps)
// Target molar mass (mtarg) must be given in g/mol
bool AngularDist::Initialize(const char* fname, double mtarg, double tgt_thickness, double beam_current){
	if(!init){
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
		
		// Need at least two points to calculate total reaction X-section
		if(num_points > 1){		
			init = true;

			com_theta = new double[num_points];
			dsigma_domega = new double[num_points];
			integral = new double[num_points];

			unsigned int index = 0;
			std::vector<double>::iterator iter1, iter2;
			for(iter1 = xvec.begin(), iter2 = yvec.begin(); iter1 != xvec.end() && iter2 != yvec.end(); iter1++, iter2++){
				com_theta[index] = *iter1;
				dsigma_domega[index] = *iter2;
				index++;
			}

			reaction_xsection = 0.0;
			integral[0] = 0.0;
			
			// Calculate the reaction cross-section from the differential cross section
			double x1, x2, y1, y2;
			for(unsigned int i = 0; i < num_points-1; i++){
				x1 = com_theta[i]*pi/180.0; y1 = dsigma_domega[i]*std::sin(x1);
				x2 = com_theta[i+1]*pi/180.0; y2 = dsigma_domega[i+1]*std::sin(x2);
				reaction_xsection += 0.5*(x2-x1)*(y2+y1);
				integral[i+1] = reaction_xsection*2*pi; // The cumulative integral
			}
			
			reaction_xsection *= 2*pi;
			//rate = avagadro*tgt_thickness*reaction_xsection*(1E-27)/(500*mtarg); // Reaction probability
			rate = 0.0;
			rate *= beam_current; // Reaction rate (pps)
			return true;
		}
		else{ 
			inFile.close();
			return false; 
		}
	}
	return false;
}

// Get a random angle from the distribution
// Returns false if a match is not found for whatever reason
bool AngularDist::Sample(double &com_angle){
	if(!init){ return false; }
	double rand_xsect = frand()*reaction_xsection;
	for(unsigned int i = 0; i < num_points-1; i++){
		if(integral[i] <= rand_xsect && rand_xsect <= integral[i+1]){ 
			com_angle = com_theta[i] + (rand_xsect-integral[i])*(com_theta[i+1]-com_theta[i])/(integral[i+1]-integral[i]);
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

// Get a random point on a circle
// spot_ is the beamspot size in m
// offset_ is the offset in the negative z-direction (in m)
// beam is a 2d vector in the xy-plane (z=0) pointing from the origin to a point inside the target beamspot
void RandomCircleUp(double spot_, double offset_, Vector3 &beam){ // Upstream of target
	double ranR = std::sqrt(frand()) * (spot_/2.0); // Random distance from the beam axis
	double ranT = 2*pi*frand(); // Random angle about the beam axis
	beam = Vector3(ranR*std::cos(ranT), ranR*std::sin(ranT), -offset_);
}

void RandomCircleDown(double spot_, double offset_, Vector3 &beam){ // Downstream of target
	double ranR = std::sqrt(frand()) * (spot_/2.0); // Random distance from the beam axis
	double ranT = 2*pi*frand(); // Random angle about the beam axis
	beam = Vector3(-ranR*std::cos(ranT), -ranR*std::sin(ranT), offset_);
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
	
	double t = 0.0, tsq = 0.0;  
	const double c0=2.515517, c1=0.802853, c2=0.010328; 
	const double d1=1.432788, d2=0.189269, d3=0.001308; 
	const double widthfact=0.424628450; 

	if(w == 0.0){ return 0.0; } 

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

/////////////////////////////////////////////////////////////////////
// Kindeux
/////////////////////////////////////////////////////////////////////
    	
// Mass values are input as AMU
// The ground state Q-value is given in MeV
// NrecoilStates_ is the number of excited states of the recoil 
// RecoilExStates_ is a pointer to an array of excitations for the recoil (in MeV)
// tgt_thickness_ is given in units of mg/cm^2
void Kindeux::Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, 
			 unsigned int NrecoilStates_, double *RecoilExStates_, double tgt_thickness_){
	if(!init){
		Mbeam = Mbeam_; Mtarg = Mtarg_; 
		Mrecoil = Mrecoil_; Meject = Meject_; 
		Qvalue = Qvalue_; tgt_thickness = tgt_thickness_; 
		NrecoilStates = NrecoilStates_;
		RecoilExStates = RecoilExStates_;
		init = true;
	}
}

// Set Kindeux to use angular distributions for calculating ejectile angles
// Returns false if attempt to load the distributions fails for any reason
bool Kindeux::SetDist(std::vector<std::string> &fnames, double total_targ_mass, double incident_beam_current){
	if(!init){ return false; }
	if(fnames.size() < NrecoilStates){
		std::cout << " Kindeux: Warning! Must have distributions for " << NrecoilStates << " excited states and the ground state\n";
		std::cout << " Kindeux: Received distributions for only " << fnames.size() << " states but expected " << NrecoilStates << std::endl;
		return false;
	}
	
	ang_dist = true;
	distributions = new AngularDist[NrecoilStates];
	Xsections = new double[NrecoilStates];

	// Load all distributions from file
	total_xsection = 0.0;
	for(unsigned int i = 0; i < NrecoilStates; i++){
		if(!distributions[i].Initialize(fnames[i].c_str(), total_targ_mass, tgt_thickness, incident_beam_current)){
			std::cout << "  Failed to load angular distribution file '" << fnames[i] << "'\n";
			ang_dist = false;
			break;
		}
		Xsections[i] = total_xsection;
		total_xsection += distributions[i].GetReactionXsection();
	}

	// Encountered some problem with one or more of the distributions
	if(!ang_dist){ 
		delete[] distributions; 
		delete[] Xsections;
		return false;
	}
	return true;
}

// Calculate recoil excitation energies
// Returns true if there was a reaction and false otherwise
bool Kindeux::get_excitations(double &recoilE, unsigned int &state){
	if(NrecoilStates == 0){
		state = 0;
		recoilE = RecoilExStates[state];
		return true;
	}
	else if(ang_dist){	
		// Angular dist weighted
		double rand_xsection = frand()*total_xsection;
		for(unsigned int i = 0; i < NrecoilStates-1; i++){
			if(rand_xsection >= Xsections[i] && rand_xsection <= Xsections[i+1]){
				// State i has been selected
				state = i;
				recoilE = RecoilExStates[state];
				return true;
			}
		}
		// rand_xsection falls in the range (Xsections[NrecoilStates-1], total_xsection]
		// State NrecoilStates-1 has been selected
		state = NrecoilStates-1;
		recoilE = RecoilExStates[state];
		return true;
	}
	else{ 
		// Isotropic
		state = (unsigned int)(frand()*NrecoilStates);
		recoilE = RecoilExStates[state];
	}
	return true;
}

// See J. B. Ball, "Kinematics II: A Non-Relativistic Kinematics Fortran Program
// to Aid Analysis of Nuclear Reaction Angular Distribution Data", ORNL-3251
bool Kindeux::FillVars(double Beam_E, double &Ejectile_E, Vector3 &Ejectile){
	unsigned int state;
	double Recoil_Ex;
	if(!get_excitations(Recoil_Ex, state)){ 
		// No reaction occured
		return false; 
	}
	
	// In the center of mass frame
	double EjectPhi, EjectTheta;
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	double temp_angle; // Ejectile angle in the center of mass frame
	if(ang_dist){
		// Sample the angular distributions for the CoM angle of the ejectile
		if(distributions[state].Sample(temp_angle)){ EjectPhi = 2*pi*frand(); } // Randomly select phi of the ejectile
		else{ UnitSphereRandom(temp_angle, EjectPhi); } // Failed to sample the distribution
	}
	else{ UnitSphereRandom(temp_angle, EjectPhi); } // Randomly select a uniformly distributed point on the unit sphere
	
	EjectTheta = std::atan2(std::sin(temp_angle),(std::cos(temp_angle)+(Vcm/VejectCoM))); // Ejectile angle in the lab
	double temp_value = std::sqrt(VejectCoM*VejectCoM-pow(Vcm*std::sin(EjectTheta),2.0));
	double Ejectile_V = Vcm*std::cos(EjectTheta); // Ejectile velocity in the lab frame
	
	if(VejectCoM >= Vcm){ 
		// Veject is single valued
		Ejectile_V += temp_value; 
	} 
	else{ 
		// Veject is double valued, so we randomly choose one of the values
		// for the velocity, and hence, the energy of the ejectile
		if(frand() >= 0.5){ Ejectile_V += temp_value; }
		else{ Ejectile_V = Ejectile_V - temp_value; }		
	}
	
	Ejectile_E = 0.5*Meject*Ejectile_V*Ejectile_V; // Ejectile energy in the lab frame
	Ejectile = Vector3(1.0, EjectTheta, EjectPhi); // Ejectile direction unit vector
	return true;
}

// Overloaded version which also calculates data for the recoil particle
bool Kindeux::FillVars(double Beam_E, double &Ejectile_E, double &Recoil_E, Vector3 &Ejectile, Vector3 &Recoil){
	unsigned int state;
	double Recoil_Ex;
	if(!get_excitations(Recoil_Ex, state)){ 
		// No reaction occured
		return false; 
	}
	
	// In the center of mass frame
	double EjectPhi, EjectTheta;
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	double temp_angle; // Ejectile and recoil angle in the center of mass frame
	if(ang_dist){
		// Sample the angular distributions for the CoM angle of the ejectile
		if(distributions[state].Sample(temp_angle)){ EjectPhi = 2*pi*frand(); } // Randomly select phi of the ejectile
		else{ UnitSphereRandom(temp_angle, EjectPhi); } // Failed to sample the distribution
	}
	else{ UnitSphereRandom(temp_angle, EjectPhi); } // Randomly select a uniformly distributed point on the unit sphere
	
	EjectTheta = std::atan2(std::sin(temp_angle),(std::cos(temp_angle)+(Vcm/VejectCoM))); // Ejectile angle in the lab
	double temp_value = std::sqrt(VejectCoM*VejectCoM-pow(Vcm*std::sin(EjectTheta),2.0));
	double Ejectile_V = Vcm*std::cos(EjectTheta); // Ejectile velocity in the lab frame
	
	if(VejectCoM >= Vcm){ 
		// Veject is single valued
		Ejectile_V += temp_value; 
	} 
	else{ 
		// Veject is double valued, so we randomly choose one of the values
		// for the velocity, and hence, the energy of the ejectile
		if(frand() >= 0.5){ Ejectile_V += temp_value; }
		else{ Ejectile_V = Ejectile_V - temp_value; }		
	}
	
	Ejectile_E = 0.5*Meject*Ejectile_V*Ejectile_V; // Ejectile energy in the lab frame
	Ejectile = Vector3(1.0, EjectTheta, EjectPhi); // Ejectile direction unit vector
	
	// Recoil calculations (now in the lab frame)
	Recoil_E = (Beam_E+Qvalue-(0.0+Recoil_Ex)) - Ejectile_E;
	Recoil = Vector3(1.0, std::asin(((std::sqrt(2*Meject*Ejectile_E))/(std::sqrt(2*Mrecoil*Recoil_E)))*std::sin(EjectTheta)), WrapValue(EjectPhi+pi,0.0,2*pi));
	return true;
}

// Convert ejectile CoM angle to Lab angle
double Kindeux::ConvertAngle2Lab(double Beam_E, double Recoil_Ex, double Eject_CoM_angle){
	// In the center of mass frame
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass	
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	
	return(std::atan2(std::sin(Eject_CoM_angle),(std::cos(Eject_CoM_angle)+(Vcm/VejectCoM)))); // Ejectile angle in the lab
}

// Print debug information
// Does nothing if angular distributions are not set
void Kindeux::Print(){
	if(ang_dist){
		for(unsigned int i = 0; i < NrecoilStates; i++){
			std::cout << "  State " << i+1 << ":";
			std::cout << " Reaction X-Section: " << distributions[i].GetReactionXsection() << " mb";
			std::cout << "\tExpected Rate: " << distributions[i].GetRate() << " pps\n";
		}
	}
}
