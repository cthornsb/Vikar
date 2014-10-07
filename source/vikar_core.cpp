// vikar_lib.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:33:21 2014

#include "../include/vikar_core.h"
#include "../include/planar.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals (for fortran commons)
/////////////////////////////////////////////////////////////////////

const double c = 2.99792458E8; // m/s
const double pi = 3.1415926540;
const double deg2rad = pi/180.0;
const double rad2deg = 180.0/pi;
const double LN2 = 0.6931471805;
const double avagadro = 6.0221413E23; // 1/mol
double eden = 0.0;
double elni = 0.0;
double avip = 0.0;
double avz = 0.0;

/////////////////////////////////////////////////////////////////////
// Vector3 Struct
/////////////////////////////////////////////////////////////////////

void Vector3::operator = (Vector3 other){
	axis[0] = other.axis[0];
	axis[1] = other.axis[1];
	axis[2] = other.axis[2];
}

// Vector addition
void Vector3::operator += (Vector3 other){
	axis[0] += other.axis[0];
	axis[1] += other.axis[1];
	axis[2] += other.axis[2];
}

// Vector subtraction
void Vector3::operator -= (Vector3 other){
	axis[0] -= other.axis[0];
	axis[1] -= other.axis[1];
	axis[2] -= other.axis[2];
}

// Scalar multiplication
void Vector3::operator *= (double scalar){
	axis[0] *= scalar;
	axis[1] *= scalar;
	axis[2] *= scalar;
}

Vector3 Vector3::operator + (Vector3 other){
	return Vector3(axis[0]+other.axis[0], axis[1]+other.axis[1], axis[2]+other.axis[2]);
}

Vector3 Vector3::operator - (Vector3 other){
	return Vector3(axis[0]-other.axis[0], axis[1]-other.axis[1], axis[2]-other.axis[2]);
}

Vector3 Vector3::operator * (double scalar){
	return Vector3(axis[0]*scalar, axis[1]*scalar, axis[2]*scalar);
}

// Dot product
double Vector3::Dot(Vector3 other){
	return (axis[0]*other.axis[0] + axis[1]*other.axis[1] + axis[2]*other.axis[2]);
}

// Cross product
Vector3 Vector3::Cross(Vector3 other){
	return Vector3((axis[1]*other.axis[2]-other.axis[1]*axis[2]),
		       (other.axis[0]*axis[2]-axis[0]*other.axis[2]),
		       (axis[0]*other.axis[1]-other.axis[0]*axis[1]));
}

// Return the length of the vector
double Vector3::Length(){
	return std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
}

double Vector3::Distance(Vector3 other){
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
	
void Vector3::Dump(){ std::cout << " " << axis[0] << ", " << axis[1] << ", " << axis[2] << std::endl; }

// Read NewVIKAR VANDLE efficiency file and load data into arrays
// Returns the number of data points which are found in the file
// Assumes the following efficiency file format for each point in file
// Lab_Energy(MeV) Efficiency
bool Efficiency::_read_eff_file(const char* fname, std::vector<double> &energy, std::vector<double> &efficiency){
	std::ifstream eff_file(fname);
	if(!eff_file.good()){ return false; }
	float values[2];

	while(true){
		for(unsigned short i = 0; i < 2; i++){ 
			eff_file >> values[i];
		}

		energy.push_back(values[0]);
		efficiency.push_back(values[1]);
		if(eff_file.eof()){ break; }
	}	
	eff_file.close();

	if(energy.size() != efficiency.size()){ return false; }
	return true;
}

Efficiency::Efficiency(){
	small_energy = NULL; small_efficiency = NULL;
	med_energy = NULL; med_efficiency = NULL;
	large_energy = NULL; large_efficiency = NULL;
	NsmallEff = 0; NmedEff = 0; NlargeEff = 0;
	init_small = false; 
	init_med = false; 
	init_large = false;
}

Efficiency::~Efficiency(){
	if(init_small){ delete[] small_energy; delete[] small_efficiency; }
	if(init_med){ delete[] med_energy; delete[] med_efficiency; }
	if(init_large){ delete[] large_energy; delete[] large_efficiency; }
}

// Load small bar efficiency data
unsigned short Efficiency::ReadSmall(const char* fname){
	std::vector<double> E, Eff;
	if(init_small || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	small_energy = new double[E.size()];
	small_efficiency = new double[E.size()];
	for(unsigned short i = 0; i < E.size(); i++){
		small_energy[i] = E[i];
		small_efficiency[i] = Eff[i];
	}
	init_small = true;

	NsmallEff = E.size();
	return NsmallEff;
}

// Load medium bar efficiency data
unsigned short Efficiency::ReadMedium(const char* fname){
	std::vector<double> E, Eff;
	if(init_med || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	med_energy = new double[E.size()];
	med_efficiency = new double[E.size()];
	for(unsigned short i = 0; i < E.size(); i++){
		med_energy[i] = E[i];
		med_efficiency[i] = Eff[i];
	}
	init_med = true;

	NmedEff = E.size();
	return NmedEff;
}

// Load medium bar efficiency data
unsigned short Efficiency::ReadLarge(const char* fname){
	std::vector<double> E, Eff;
	if(init_large || !_read_eff_file(fname, E, Eff)){ return 0; }
	
	// Generate the efficiency arrays
	large_energy = new double[E.size()];
	large_efficiency = new double[E.size()];
	for(unsigned short i = 0; i < E.size(); i++){
		large_energy[i] = E[i];
		large_efficiency[i] = Eff[i];
	}
	init_large = true;

	NlargeEff = E.size();
	return NlargeEff;
}

// Return the interpolated value for the input energy
double Efficiency::GetSmallEfficiency(double Energy){
	if(!init_small || NsmallEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < small_energy[0]){ efficiency = small_efficiency[0]; }
	else if(Energy > small_energy[NsmallEff-1]){ efficiency = small_efficiency[NsmallEff-1]; }
	else{
		for(unsigned short i = 1; i < NsmallEff; i++){
			if(Energy >= small_energy[i-1] && Energy <= small_energy[i]){
				efficiency = Interpolate(small_energy[i-1],small_efficiency[i-1],
							 small_energy[i],small_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
}

// Return the interpolated value for the input energy
double Efficiency::GetMediumEfficiency(double Energy){
	if(!init_med || NmedEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < med_energy[0]){ efficiency = med_efficiency[0]; }
	else if(Energy > med_energy[NmedEff-1]){ efficiency = med_efficiency[NmedEff-1]; }
	else{
		for(unsigned short i = 1; i < NmedEff; i++){
			if(Energy >= med_energy[i-1] && Energy <= med_energy[i]){
				efficiency = Interpolate(med_energy[i-1],med_efficiency[i-1],
							 med_energy[i],med_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
}	

// Return the interpolated value for the input energy
double Efficiency::GetLargeEfficiency(double Energy){
	if(!init_large || NlargeEff == 0){ return 1.0; }
	
	double efficiency = 0.0;
	if(Energy < large_energy[0]){ efficiency = large_efficiency[0]; }
	else if(Energy > large_energy[NlargeEff-1]){ efficiency = large_efficiency[NlargeEff-1]; }
	else{
		for(unsigned short i = 1; i < NlargeEff; i++){
			if(Energy >= large_energy[i-1] && Energy <= large_energy[i]){
				efficiency = Interpolate(large_energy[i-1],large_efficiency[i-1],
							 large_energy[i],large_efficiency[i],Energy);
				break;
			}
		}
	}
	return efficiency;
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
		while(true){
			if(inFile.eof()){ break; }
			inFile >> x >> y;
			com_theta.push_back(x*deg2rad);
			dsigma_domega.push_back(y);
			num_points++;
		}
		
		// Need at least two points to calculate total reaction X-section
		if(num_points > 1){		
			init = true;
			reaction_xsection = 0.0;
			integral.push_back(0.0);
			
			// Calculate the reaction cross-section from the differential cross section
			for(unsigned int i = 0; i < num_points-1; i++){
				//reaction_xsection += (com_theta[i+1]-com_theta[i])*(dsigma_domega[i]+dsigma_domega[i+1])/2.0; // Trapezoidal integration
				reaction_xsection += (com_theta[i+1]-com_theta[i])*(dsigma_domega[i]*std::sin(com_theta[i])+dsigma_domega[i+1]*std::sin(com_theta[i+1]))/2.0; // Trapezoidal integration
				integral.push_back(reaction_xsection*2*pi);
			}
			reaction_xsection *= 2*pi;
			rate = avagadro*tgt_thickness*reaction_xsection*(1E-27)/(500*mtarg); // Reaction probability
			rate *= beam_current; // Reaction rate (pps)
		}
		else{ 
			inFile.close();
			return false; 
		}
	}
	return true;
}

// Get a random angle from the distribution
// Returns 0.0 if a match is not found for whatever reason
double AngularDist::Sample(){
	double rand_xsect = frand()*reaction_xsection;
	for(unsigned int i = 0; i < num_points-1; i++){
		if(integral[i] <= rand_xsect && rand_xsect <= integral[i+1]){ return (com_theta[i]+com_theta[i+1])/2.0; }
	}
	return 0.0;
}

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

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

// Convert char* to short
short atos(const char* input){
	return short(atol(input));
}

// Sample a point on the surface of the unit sphere
void UnitRandom(Vector3 &vec){
	double u = 2*frand()-1;
	double theta = 2*pi*frand();
	vec.axis[0] = std::sqrt(1-u*u)*std::cos(theta);
	vec.axis[1] = std::sqrt(1-u*u)*std::sin(theta);
	vec.axis[2] = u;
}

void UnitRandom(double &theta, double &phi){
	phi = 2*pi*frand();
	theta = std::acos(2*frand()-1);
}

// Calculate proper bar spacing for a wall of VANDLE bars
// Leave half gaps at either edge for clearance to other walls
double BarSpacing(double total_width, double bar_width, unsigned short num_bars){
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
unsigned short GetLines(const char* input){
	std::ifstream file(input); 
	unsigned short count = 0;
	std::string line;
	while(std::getline(file,line)){ 
		count++; 
	}
	file.close();
	
	return count;
}

// Convert arbitrary input to string
std::string to_str(double input){
	std::stringstream output;
	output << input;
	return output.str();
}

// Linearly interpolate between points
// Return the value y = f(x)
double Interpolate(double x1, double y1, double x2, double y2, double x){
	return ((y2-y1)/(x2-x1))*(x-x1)+y1;
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

void Cart2Sphere(Vector3 cart, Vector3 &sphere){
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
	
	/*x = r*std::sin(theta)*std::cos(phi); 
	y = r*std::sin(theta)*std::sin(phi); 
	z = r*std::cos(theta);*/
	x = r*std::sin(theta)*std::cos(phi); 
	y = r*std::sin(theta)*std::sin(phi); 
	z = r*std::cos(theta);
}

void Sphere2Cart(double r, double theta, double phi, Vector3 &cart){
	cart.axis[0] = r*std::sin(theta)*std::cos(phi);
	cart.axis[1] = r*std::sin(theta)*std::sin(phi); 
	cart.axis[2] = r*std::cos(theta);
}

void Sphere2Cart(Vector3 sphere, Vector3 &cart){
	cart.axis[0] = sphere.axis[0]*std::sin(sphere.axis[1])*std::cos(sphere.axis[2]);
	cart.axis[1] = sphere.axis[0]*std::sin(sphere.axis[1])*std::sin(sphere.axis[2]); 
	cart.axis[2] = sphere.axis[0]*std::cos(sphere.axis[1]);
}

/////////////////////////////////////////////////////////////////////
// dedx.f
/////////////////////////////////////////////////////////////////////

double beta2(double e, double em){
	// to calculate beta**2 where beta is the speed of a particle
	// relative to the speed of light
	//
	// input - e energy in mev
	// em rest mass in mev
	
	//double r = e/em + 1.0; 
	return 1.0-(1.0/(pow(e/em + 1.0, 2.0)));
} 

double btoep(double dbsq){
	// to calculate the energy of a proton given beta**2, where
	// beta is the speed of the proton relative to the speed of light
	// input:-   dbsq - beta**2 (double precision)
	 
	//const double emp = 938.25920; 
	//double d = std::sqrt(1.0/(1.0-dbsq))-1.0; 
	return 938.25920*(std::sqrt(1.0/(1.0-dbsq))-1.0);
} 

double zeff(double z, double beta){
	//          to calculate the effective charge of a particle in
	//     stopping power calculations(spar-armstrong & chandler-ornl-
	//     4869 [1973])
	//
	//     input - z nominal charge of nucleus
	//             beta speed of the particle relative to the speed of light
	 
	double zp = pow(z, 0.6666670); 
	if(beta <= (0.070*zp)){ return z*(1.00 - std::exp(-125.00*beta/zp)); }
	else{ return z; }
} 

double deff(double e){ 
	//          to calculate the density-effect correction term in stopping
	//     power calculations (spar-armstrong & chandler-ornl-4869 [1973])
	//
	//     input - e energy in mev
	//             elni astd::log(mean ionization potential of medium (mev) )
	//             eden electron density of stopping medium
	 
	
	double del; 
	//const double emp1 = 0.00106580356; 
	//const double emp = 938.25920;
	 
	del = e*0.00106580356; 
	del = std::log(del) + std::log(del+2.00); 
	del += std::log(1.378E-9*eden) - 2.00*elni - 1.00; 
	
	if(del <= 0.00){ return 0.00; }
	else{ return del; }
}

double shell(double e){
	//          to calculate the shell correction term in stopping power
	//     calculations(spar-armstrong & chandler-ornl-4869 [1973])
	//
	//     input - e energy in mev
	//             avip mean ionization potential of stopping medium (mev)
	//             avz mean atomic number of stopping medium
	 
	double be2; 
	double gnu2; 
	double f1, f2, xl, xl1, xl2; 
	double output, x, xlog;
	
	const double a1 = 0.422377E-6, a2 = 3.858019E-9, b1 = 0.304043E-7; 
	const double b2 = -0.1667989E-9, c1 = -0.38106E-9, c2 = 0.157955E-11; 
	const double emp = 938.25920, zal = 12.950, zh2o = 3.340; 
	const double emp1 = 0.00106580356; 
	
	const double p1 = 4.774248E-4, p2 = 1.143478E-3, p3 = -5.63392E-2; 
	const double p4 = 4.763953E-1, p5 = 4.844536E-1; 
	const double w1 = -1.819954E-6, w2 = -2.232760E-5, w3 = 1.219912E-4; 
	const double w4 = 1.837873E-3, w5 = -4.457574E-3, w6 = -6.837103E-2; 
	const double w7 = 5.266586E-1, w8 = 3.743715E-1; 
	const double const1 = 0.021769515; 

	// originally (E >= 8.0) SDP
	if(e >= 2008.0){
		gnu2 = 1.0/((e*emp1)*(e*emp1+2.00)); 
		f1 =  gnu2*(a1+gnu2*(b1+c1*gnu2)); 
		f2 =  gnu2*(a2+gnu2*(b2+c2*gnu2)); 
		return avip*avip*(f1 + f2*avip)/avz; 
	} 

	be2 = beta2(e, emp); 
	output = const1 + std::log(be2) - std::log(avip); 
	x = 18769.00*be2/avz; 
	xlog=std::log(x);
	xl1 = 0.0;
	
	if(avz > zal){
		xl = p5 + xlog*(p4+xlog*(p3+xlog*(p2+xlog*p1)));
		xl = std::exp(xl);
		if(avz > zal){ output = output - xl; }
		else{
			xl2 = xl;
			xl = xl1 + (avz-zh2o)/(zal-zh2o)*(xl2-xl1);
			output = output - xl;
		}
	}
	else{
		xl1 = w8 + xlog*(w7+xlog*(w6+xlog*(w5+xlog*(w4+xlog*(w3+xlog*(w2+xlog*w1))))));
		xl1 = std::exp(xl1);
		if(avz > zh2o){
			xl = p5 + xlog*(p4+xlog*(p3+xlog*(p2+xlog*p1)));
			xl = std::exp(xl);
			if(avz > zal){ output = output - xl; }
			else{
				xl2 = xl;
				xl = xl1 + (avz-zh2o)/(zal-zh2o)*(xl2-xl1);
				output = output - xl;
			}	
		}
		else{
			xl = xl1;
			output = output - xl;
		}
	} 
	
	return output;
} 

double dedxp(double enrgy, double db2, double beta){
	// only valid for beta > 0.0046 : 10 kev protons
	
	//const double emass = 938.2592;
	return eden*0.5099147*(pow(zeff(1.0,beta), 2.0)/db2)*((log(1.022008*db2/(1.0-db2))-db2)-elni-shell(enrgy)-deff(enrgy)*0.5);
	//double ze = zeff(1.0,beta);
	//dsp *= (ze*ze)/db2;
	//dsp *= pow(zeff(1.0,beta), 2.0)/db2;
	//double d = log(1.022008*db2/(1.0-db2))-db2;
	//double delta = deff(enrgy);
	//double coz = shell(enrgy);
	//dsp *= (log(1.022008*db2/(1.0-db2))-db2)-elni-shell(enrgy)-deff(enrgy)*0.5;
	//return dsp;
} 

double dedx(double emass, double epart, double zpart){ 
	// only valid for beta > 0.0046*z**(1/3)
	  
	double db2 = beta2(epart,emass);
	double beta = std::sqrt(db2);
	//double pe = btoep(db2);
	//double gog = zeff(zpart,beta)/zeff(1.0,beta);
	return pow(zeff(zpart,beta)/zeff(1.0,beta), 2.0)*dedxp(btoep(db2),db2,beta);
} 

double range(double dx, double emass, double epart, double zpart, double sp){
	double spe = sp; 
	double ed = epart; 
	double es = epart; 
	double dxt = dx; 
	double dxs = dxt; 
	double f = 1.00; 
	double g = 1.00; 
	double er2 = 1.00; 
	double e, r, er1, jf;
	e = 0.0; r = 0.0; er1 = 0.0; jf = 0.0;

	int count = 0;
	while(count < 10000){
		ed = ed-de(dxt,emass,es,zpart,spe);
		e = ed;
		if(e <= 0.0000001){ break; }
		
		r = r+dxt;
		es = e;
		spe = dedx(emass,es,zpart);
		er1 = epart/e;
		
		if(std::abs(er1/er2 - 1.0) < 0.01505){
			g *= 1.10;
			dxt = g*dxs;
			dxs = dxt;
		}
		else{
			jf = min(er1, 32.0);
			f = 1.0/jf;
		}
		
		er2 = er1;
		dxt = f*dxs;
		count++;
	} 

	//if(count >= 10000){ std::cout << " Warning: Value did not converge!\n"; }
	if(count >= 10000){ 
		std::cout << ed << " " << dxt << " " << emass << " " << es << " " << zpart << " " << spe << " " << e << " " << r << std::endl; 
	}
	return r+es/spe;
} 

double range2(double dx, double emass, double epart, double zpart, double sp){
	double spe = sp;
	double ed = epart;
	double r = 0.0;
	double es = epart;
	double dxt = dx;
	double dxs = dxt;
	double f = 1.0;
	double g = 1.0;
	double er2 = 1.0;
	double e, er1, jf;
	e = 0.0; er1 = 0.0; jf = 0.0;
      
    top:
	ed = ed - de(dxt,emass,es,zpart,spe);
	e = ed;
	if(e <= 0.0000001){ goto bottom; }
	r = r + dxt;
	es = e;
	spe = dedx(emass,es,zpart);
	er1 = epart/e;

	if(abs(er1/er2 - 1.0) < 0.01505){
		g = g*1.10;
		dxt = g*dxs;
		dxs = dxt;
	}
	else{
		jf = min(er1,32);
		f = 1.0/jf;
	}
	
	er2 = er1;
	dxt = f*dxs;
	std::cout << ed << " " << dxt << " " << emass << " " << es << " " << zpart << " " << spe << " " << e << " " << r << std::endl;
	goto top;
	
    bottom:
	return r+es/spe;
}

/////////////////////////////////////////////////////////////////////
// ncdedx.f
/////////////////////////////////////////////////////////////////////

double algip(double z){ 
	//          to calculate alog(ionization potential) for an element
	//     of atomic number z
	//
	//      n.b. ionization potl in mev !!!!!!
	double pot[13] = {18.7,42.0,39.0,60.0,68.0,78.0,99.5,98.5,117.0,140.0,150.0,157.0,163.0};
	double iz, potl;
	
	iz = z + 0.050;
	if(iz > 12){ potl = 9.760*z + 58.80/(pow(z, 0.190)); }
	else{ potl = pot[short(iz)-1]; }
	return std::log(potl * 1.0e-6);
} 

void ncdedx(double tgtdens, double atarget, double ztarget, double abeam, double zbeam, double energy, 
	    double &dedxmg, double &tgtionpot, double &rangemg){ 
	double mchem[10], zmed[10], amua[10];
	double nmed = 1.0;
	double avden = tgtdens;
	amua[0] = atarget;
	zmed[0] = ztarget;
	mchem[0] = 1.0;
	double amum = mchem[0]*amua[0];
	double denm = avden/(amum*1.660543);
	
	double sumn = 0.0;
	double sumnz = 0.0;
	double sumnzi = 0.0;
	double en,enz,enzi;
	for(unsigned short i = 0; i < nmed; i++){
		en = mchem[i]*denm;
		enz = en*zmed[i];
		enzi = enz*algip(zmed[i]);
		sumn += en;
		sumnz += enz;
		sumnzi += enzi;
	}
	
	elni = sumnzi/sumnz; //global
	avip = std::exp(elni); //global
	tgtionpot = avip*1000000.0; //return (good)
	eden = sumnz; //global
	avz = sumnz/sumn; //global
	double emass = abeam*931.4812;
	double epart = energy;
	double dx = 0.5/dedx(emass,epart,zbeam);
	double sp = dedx(emass,epart,zbeam);
	double r = range2(dx,emass,epart,zbeam,sp);
	dedxmg = (sp*0.001)/avden; //return (good)
	rangemg = r*avden*1000.0; //return (off)
} 

/////////////////////////////////////////////////////////////////////
// de.f
/////////////////////////////////////////////////////////////////////

double de(double dx, double emass, double epart, double zpart, double sp){ 
	// to calculate energy lost over dx cms assuming a quadratic relationship between e & x.
	double deltae = sp*dx; 
	double enew = epart-deltae; 
	//std::cout << " de: " << sp << " " << dx << " " << epart << " " << deltae << std::endl;
	
	if(enew <= 0.0){ return epart; }
	
	double spges = dedx(emass, enew, zpart); 
	return dx*(0.750*sp+(0.250*spges/sp)*spges); 
} 

/////////////////////////////////////////////////////////////////////
// linear.f
/////////////////////////////////////////////////////////////////////

double linear(double xmin, double xmax, double ymin, double ymax, double x){
	// linear 1.0 written by S.D.Pain on 24/11/2004
	// Function for linear interpolation between two points

	double grad = (ymax-ymin)/(xmax-xmin); 
	double cint = ymin - grad*xmin; 
	return x*grad + cint; 
} 

/////////////////////////////////////////////////////////////////////
// momentum.f
/////////////////////////////////////////////////////////////////////

double momentum(double energy, double mass){
	// momentum 1.0 written by S.D.Pain on 27/01/2005
	// Function to calculate the momentum of a body with 'energy' and 'mass'

	return sqrt(2.0*energy*mass); 
} 

/////////////////////////////////////////////////////////////////////
// radlength.f
/////////////////////////////////////////////////////////////////////

double radlength(unsigned short A, unsigned short Z){
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

// rndgauss1 : Generate random numbers with a Gaussian distribution *
void rndgauss1(double &u, double &x, double &f, double &c, double &s){ 
	//     a subroutine to calculate one random deviate
	//     for a normal distribution with centroid c and
	//     standard deviation s
	//
	//     see chapter 26 of Abramowitz and Stegun
	//      equn 26.2.23
	//
	//      on exit u contains the random number ( deviate)
	//       -calculated by call to drand
	//
	//      on exit  x contains the deviate with a normal distribution
	//
	//      on exit the  f contains the probability function ( a gaussian)
	//
	//      of the form f(x) = 1/(s *sqrt(2*pi)) * exp(- (x-c)**2/(2*s**2))
	//
	//      The physical problem is  f(x) = u where u is a random deviate
	//      and we wish to find  x  = f**-1 ( u)
	//
	//      an inverse Chebyshev polynomial expansion is used
	//
	//      ***************************************************************
	//
	//      initial calculations assuming c=0.0, s=1.0
	//
	//       get random probability , u, in the range 0 < u <= 0.5
	//
	//     u=drand(0) *0.50
	//     usq=u*u
	//     if (usq.lt.1.0d-60) usq=1.0d-60
	//     t= sqrt(log(1.00/usq))
	//     tsq=t*t
	//     tcube=tsq*t
	
	const double c0=2.5155170, c1=0.8028530, c2=0.0103280; 
	const double d1=1.4327880, d2=0.1892690, d3=0.0013080; 
	const double rcpsqr2pi=0.398942280, sqrt2=1.4142135620; 

	double t, tsq; 

	u = frand(); 
	if(u > 0.5){ u = u-0.5; }
	if(u < 1e-30){
		t = 11.7539400024; 
		tsq = 138.1550558; 
	} 
	else{ 
		tsq = -std::log(u*u); 
		t = std::sqrt(tsq); 
	} 

	//       compute inverse by equn 26.2.23
	x = t-(c0+c1*t+c2*tsq)/(1.00+d1*t+(d2+d3*t)*tsq); 

	//     now randomize x in positive and negative direction
	//     x=x* (2* nshort(drand(0)) -1)
	if(u > 0.5){ x=-x; }

	//     compute function
	f = rcpsqr2pi*std::exp(-(x*x)); 

	//     now correct for centroid and standard deviation
	x = sqrt2*s*x+c; 
	f = f/s; 
} 

/////////////////////////////////////////////////////////////////////
// velocity.f
/////////////////////////////////////////////////////////////////////

double velocity(double energy, double mass){ 
	// strag_targ 1.0 written by S.D.Pain on 20/11/2004
	// Function to calculate the velocity of a body with 'energy' and 'mass'

	return std::sqrt(2.0*energy/mass); 
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

	theta = 13.6/(velocity(energy,A)*momentum(energy,A))*Z*std::sqrt(thickness/X)*(1.0+0.038*std::log(thickness/X)); 
	theta = theta*std::sqrt(2.0); 
} 

/////////////////////////////////////////////////////////////////////
// transform.f
/////////////////////////////////////////////////////////////////////

void transform(double theta1, double phi1, double theta2, double phi2, double &theta, double &phi){ 
	// transform 2.0 written by S.D.Pain on 4/03/2005
	//
	// Subroutine for transforming the a spherical polar vector
	// from one refernce frame to another.
	// (theta2,phi2) is a vector in the master frame
	// (theta1,phi1) is measured relative to (theta2,phi2).
	// (theta,phi) is (theta1,phi1) in the master frame

	double term1, term2, temp, x1, y1, x2, y2, x, y; 
	double beamX, beamY, beamZ, dummy, pi; 
	double dumtheta1, dumphi1, dumtheta2, dumphi2; 
	bool swap; 
	
	swap = false; 
	
	pi = 3.1415926540; 
	dummy = 1.0; 
	
	// copy the input angles to different variables, and use the copies in
	// the subroutine, as they get modified.
	dumtheta1 = theta1; 
	dumphi1 = phi1; 
	dumtheta2 = theta2; 
	dumphi2 = phi2; 
	
	// Check whether the vector is pointing backward of 90 degrees (polar)
	// If so, reflect its direction around, so that it points forwards.
	// The transformation can then be computed, and the vector reflected
	// back again. This avoids edge-of-the-world effects.
	
	if (dumtheta1 > (0.5*pi)){
		swap = true; 
		
		// Doesn't appear that any transformation is needed here - perhaps
		// worth checking, though...
		//        call Sphere2Cart(dummy,theta2,phi2,beamX,beamY,beamZ)
		//        beamX = -beamX
		//        beamY = -beamY
		//        beamZ = -beamZ
		//        call Cart2Sphere(beamX,beamY,beamZ,dummy,theta2,phi2)
		
		Sphere2Cart(dummy,dumtheta1,dumphi1,beamX,beamY,beamZ); 
		beamX = -beamX; 
		beamY = -beamY; 
		beamZ = -beamZ; 
		Cart2Sphere(beamX,beamY,beamZ,dummy,dumtheta1,dumphi1); 
	} 
	
	// Calculate the total angle between the two vectors. This is the
	// effective polar angle.
	term1 = dumtheta1 + dumtheta2*(cos(dumphi2-dumphi1)); 
	term2 = dumtheta2*(sin((dumphi2-dumphi1))); 
	theta = sqrt(pow(term1, 2 )+pow( term2, 2) ); 
	
	x1 = dumtheta1*sin(dumphi1); 
	y1 = dumtheta1*cos(dumphi1); 
	x2 = dumtheta2*sin(dumphi2); 
	y2 = dumtheta2*cos(dumphi2); 
	
	x = x1+x2; 
	y = y1+y2; 
	
	temp = x/(sqrt(pow(x, 2)+pow(y, 2)) ); 
	phi = asin(temp); 
	
	if (x>=0.0){ phi = acos(y/(std::sqrt(pow(x, 2)+pow(y, 2)))); } 
	else{ phi = 2.0*3.14159-acos(y/(sqrt(pow(x, 2)+pow(y, 2)))); }
	
	// Keeps theta & phi within limits. Not needed, with reflection
	// procedure
	//      if(theta.gt.3.14159)then
	//       theta = 2.0*3.14159-theta
	//        phi = phi+3.14159
	//        if(phi.gt.2.0*3.14159) phi = phi-2.0*3.14159
	//      endif
	
	// If a reflection was made, reflect back again.
	if (swap){
		Sphere2Cart(dummy,theta,phi,beamX,beamY,beamZ); 
		beamX = -beamX; 
		beamY = -beamY; 
		beamZ = -beamZ; 
		Cart2Sphere(beamX,beamY,beamZ,dummy,theta,phi); 
	} 
} 

/////////////////////////////////////////////////////////////////////
// strag_targ.f
/////////////////////////////////////////////////////////////////////

void strag_targ(double A, double Z, double targ_thick, double theta_old, double phi_old, double energy, double &theta_new, double &phi_new, double X){ 
	// strag_targ 1.0 written by S.D.Pain on 20/01/2004
	//
	// strag_targ 1.1 modified by S.D.Pain on 7/03/2005
	// to untilise transform2.0
	//
	// Subroutine to calculate the angular straggling of an ion in
	// the target. The A,Z of the ion are read in, along with the
	// theta,phi and energy of the ion. The average radiation length
	// of the target material is read in as X.
	//
	// The calculation of the width of the scattering distribution
	// is calculated by straggleA. A Gaussian weighted scattering angle
	// based on this width is calculated, using rndgauss0
	//
	// The new theta,phi to which the ion is scattered is returned.
	
	double theta_scatW; 
	double theta_scat, phi_scat; 
	
	// Calculate the straggling width
	straggleA(theta_scatW, energy, Z, A, targ_thick, X); 
	
	// Select the scattering angle of the ion wrt its initial direction
	theta_scat = rndgauss0(theta_scatW); 
	theta_scat = std::sqrt(pow(theta_scat, 2)*2.0); 
	phi_scat = frand()*2.0*pi; 
	
	// Determine the std::absolute lab angle to which the ion is scattered
	transform(theta_old, phi_old, theta_scat, phi_scat, theta_new, phi_new); 
}

/////////////////////////////////////////////////////////////////////
// targ_thick.f
/////////////////////////////////////////////////////////////////////

void targ_thick(double theta_in, double phi_in, double thicknessZ, double depth, double theta_targ, double &thickness){ 
	// targ_thick 1.0 written by S.D.Pain on 2/2/2005
	//
	//   theta_in = polar angle of particle in lab
	//   phi_in = azimuthal angle of particle in lab
	//   thicknessZ = total target thickness (mg/cm^2)
	//   depth = depth through target for interaction, perpendicular to target plane
	//   theta_targ = angle of target (rotated clockwise about y axis, viewed from above)
	//   thickness = thickness of material the particle must pass through

	double x, y, z; 
	double dummy, pi; 
	double newx, newy, newz; 
	double dummyr, dummytheta, dummyphi; 
	pi = 3.14159; 
	
	// Convert the ion's direction to cartesian coordinates
	dummy = 1.0; 
	Sphere2Cart(dummy,theta_in,phi_in,x,y,z); 
	//length = sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2)); 
	
	// Rotate the ion's vector so it is measured wrt to the target
	newx = x*cos(theta_targ) - z*sin(theta_targ); 
	newz = z*cos(theta_targ) + x*sin(theta_targ); 
	newy = y; 
	
	// Convert the ion's vector back to spherical polars, now wrt the
	// target
	Cart2Sphere(newx, newy, newz, dummyr, dummytheta, dummyphi); 
	
	// Calculate the thickness seen by the ion. Bloody marvelous.
	if(dummytheta <= (0.5*pi)){ thickness = (thicknessZ-depth)/std::cos(dummytheta); }
	else{ thickness = (-depth)/std::cos(dummytheta); }
	
	// This line was added to account for an occasion where dummytheta was greater
	// than 0.5*pi (only just), but for some reason, cos(dummytheta) was positive
	// which gave a negative thickness.
	if(thickness < 0){ thickness = -thickness; }
}

/////////////////////////////////////////////////////////////////////
// unitV.f
/////////////////////////////////////////////////////////////////////

void unitV(double xl, double yl, double zl, double &x, double &y, double &z, double &length){ 
	// unitV 1.0 written by S.D.Pain on 20/11/2004
	//
	// Subroutine to read in a cartesian vector (xl,yl,zl)
	// and calculate and return its unit vector (x,y,z) and length
	
	// Put in check of length > 0    ?????
	length = std::sqrt(pow(xl, 2)+pow(yl, 2)+pow(zl, 2)); 
	x = xl/length; 
	y = yl/length; 
	z = zl/length; 
}

/////////////////////////////////////////////////////////////////////
// AngDist_read.f
/////////////////////////////////////////////////////////////////////

void AngDist_read(std::string fName, unsigned short &Npoints, double *angle, double *integral, double &max_integral){ 
	// AngDist_read 1.0 written by S.D.Pain on 5/05/2006
	//
	// Subroutine for reading in an angular distribution profile from
	// from fName.
	// An angular distribution profile is a cumulative integration of
	// the angular distribution.
	// The input file should be of the form:
	// [angle (deg)] [cumulative integral(0-angle)]
	// where the angles must span the range 0 to 180 degrees.	
	
	unsigned short i; 
	max_integral = 0.0; 
	// DetSet_stat stores error status - T = good, F = bad

	std::ifstream file10(fName.c_str()); 
	Npoints = 0; // Zero the points counter
	
	// Read in the main data points from the SRIM output file
	while(!file10.eof()){
		Npoints = Npoints+1; 
		file10 >> angle[Npoints] >> integral[Npoints]; 
	} // Read in data loop
	
	if(Npoints > 0){ Npoints = Npoints-1; }
	file10.close(); 
	
	for (i = 0; i < Npoints; i++){
		max_integral = integral[Npoints]; 
	} 
}

/////////////////////////////////////////////////////////////////////
// srim-read.f
/////////////////////////////////////////////////////////////////////

void SRIMread(std::string fName, bool &SRIM_stat, unsigned short &Npoints, double *energy, double *dedx, double *range, double *longitude, double *latitude, bool convert){ 
	// SRIMread 1.0 written by S.D.Pain on 24/11/2004
	//
	// Subroutine for reading in data from a SRIM output file fName
	// The material density is read in, and all data in the main table
	// are read in.
	// Conversions are made to ensure all measurements are in um
	// Ranges (and widths) are converted from um to mg/cm^2,
	// using the density read from the SRIM file
	// SRIM_stat provides limited error reporting (!)
	// SRIM_stat onwards are retunred variables	
 
	std::string dummyC, unit_E, unit_range, unit_longitude, unit_latitude; 
	double dedxE, dedxN, density, conv;
	std::string junk;
	
	// SRIM_stat stores error status - T = good, F = bad
	SRIM_stat = true; 
	
	std::ifstream file10(fName.c_str()); 
	for(short i = 1; i <= 10; i++){
		file10 >> junk; 
	} 
	
	// Read in density in g/cm^3
	file10 >> dummyC >> dummyC >> dummyC >> density; 
	density = density*1000.0; // *16.0/14.0 // Convert density to mg/cm^3
	conv = 1.0e-4*density; 
	for(short i = 1; i <= 4; i++){
		file10 >> junk; 
	} 
	
	while(dummyC == "------"){
		file10 >> dummyC; 
	} 
	
	Npoints = 0; // Zero the data pounsigned short counter
	
	// Read in the main data points from the SRIM output file
	while(!file10.eof() && SRIM_stat){
		Npoints=Npoints+1; 
		file10 >> energy[Npoints] >> unit_E >> dedxE >> dedxN >> range[Npoints] >> unit_range;
		file10 >> longitude[Npoints] >> unit_longitude >> latitude[Npoints] >> unit_latitude; 
		
		// Add the dedx for electric and nuclear effects
		dedx[Npoints] = dedxE + dedxN; 
		
		// Make sure the energies are in MeV
		if(unit_E == "eV"){ energy[Npoints] = energy[Npoints]/1000000.0; }
		else if(unit_E == "keV"){ energy[Npoints] = energy[Npoints]/1000.0; }
		else if(unit_E == "MeV"){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the range values are in um
		if(unit_range ==  "A"){ range[Npoints] = range[Npoints]*0.0001; }
		else if(unit_range ==  "mm"){ range[Npoints] = range[Npoints]*1000.0; }
		else if(unit_range ==  "um"){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the longitude values are in um
		if(unit_longitude ==  "A"){ longitude[Npoints] = longitude[Npoints]*0.0001; }
		else if(unit_longitude ==  "mm"){ longitude[Npoints] = longitude[Npoints]*1000.0; }
		else if(unit_longitude ==  "um"){  } 
		else{ SRIM_stat = false; } 
		
		// Make sure the latitude values are in um
		if(unit_latitude ==  "A"){ latitude[Npoints] = latitude[Npoints]*0.0001; }
		else if(unit_latitude ==  "mm"){ latitude[Npoints] = latitude[Npoints]*1000.0; }
		else if(unit_latitude ==  "um"){  } 
		else{ SRIM_stat = false; } 
		
		// Convert from length to mg/cm^2 if necessary
		if (convert){
			range[Npoints] = range[Npoints]*conv; 
			longitude[Npoints] = longitude[Npoints]*conv; 
			latitude[Npoints] = latitude[Npoints]*conv; 
		} 
	} // Read in data loop
	
	Npoints = Npoints-1; 
	file10.close(); 

	if(!SRIM_stat){
		std::cout << "Arse Biscuits!!";
	} 
} 

/////////////////////////////////////////////////////////////////////
// kindeux.f
/////////////////////////////////////////////////////////////////////
    	
// Mass values are input as AMU
// The ground state Q-value is given in MeV
// NrecoilStates_ is the number of excited states of the recoil 
// RecoilExStates_ is a pointer to an array of excitations for the recoil (in MeV)
// tgt_thickness_ is given in units of mg/cm^2
void Kindeux::Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_, 
			 unsigned short NrecoilStates_, double *RecoilExStates_, double tgt_thickness_){
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
	if(init){
		if(fnames.size() < NrecoilStates){
			std::cout << " Kindeux: Warning! Must have distributions for " << NrecoilStates << " excited states and the ground state\n";
			std::cout << " Kindeux: Received distributions for only " << fnames.size() << " states but expected " << NrecoilStates << std::endl;
			return false;
		}
		
		NDist = fnames.size();
		ang_dist = true;
		distributions = new AngularDist[NDist];
	
		// Load all distributions from file
		for(unsigned short i = 0; i < NDist; i++){
			if(!distributions[i].Initialize(fnames[i].c_str(), total_targ_mass, tgt_thickness, incident_beam_current)){
				ang_dist = false;
				break;
			}
		}
	
		// Encountered some problem with one or more of the distributions
		if(!ang_dist){ 
			delete[] distributions; 
			return false;
		}
		return true;
	}
	return false;
}

// Calculate recoil excitation energies
// Returns true if there was a reaction and false otherwise
bool Kindeux::get_excitations(double& recoil){
	if(NrecoilStates == 0){
		recoil = RecoilExStates[0];
		return true;
	}
	else if(ang_dist){	
		// Angular dist weighted
		// This allows for the possibility that no reaction occurs
		// Because of this, we can actually calculate the reaction rate
		recoil = RecoilExStates[int(frand()*NrecoilStates)]; // Temporary, for testing
		if(frand() >= 0.5){ return true; }
		else{ return false; }
	}
	else{ 
		// Isotropic
		recoil = RecoilExStates[int(frand()*NrecoilStates)];
		return true;
	}
}

// See J. B. Ball, "Kinematics II: A Non-Relativistic Kinematics Fortran Program
// to Aid Analysis of Nuclear Reaction Angular Distribution Data", ORNL-3251
bool Kindeux::FillVars(double Beam_E, double theta_beam, double phi_beam, double &Ejectile_E, Vector3 &Ejectile){
	double Recoil_Ex;
	if(!get_excitations(Recoil_Ex)){ 
		// No reaction occured
		return false; 
	}
	
	// In the center of mass frame
	double EjectPhi, EjectTheta;
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	double temp_angle; // Ejectile angle in the center of mass frame
	UnitRandom(temp_angle, EjectPhi); // Randomly select a uniformly distributed point on the unit sphere
	
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
bool Kindeux::FillVars(double Beam_E, double theta_beam, double phi_beam, double &Ejectile_E, double &Recoil_E, Vector3 &Ejectile, Vector3 &Recoil){
	double Recoil_Ex;
	if(!get_excitations(Recoil_Ex)){ 
		// No reaction occured
		return false; 
	}
	
	// In the center of mass frame
	double EjectPhi, EjectTheta;
	double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass
	double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
	
	double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+Recoil_Ex))); // Ejectile CoM velocity after reaction
	double temp_angle; // Ejectile and recoil angle in the center of mass frame
	UnitRandom(temp_angle, EjectPhi); // Randomly select a uniformly distributed point on the unit sphere
	
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

// Overloaded version which also calculates data for the recoil particle
bool Kindeux::FillVars(double Beam_E, double beam_theta_com, double *EjectTheta, double *RecoilTheta, double *RecoilTheta2, double *Ejectile_E, double *Ejectile_E2, double *Recoil_E, double *Recoil_E2){
	for(unsigned int i = 0; i < NrecoilStates; i++){
		// In the center of mass frame
		double Vcm = std::sqrt(2.0*Mbeam*Beam_E)/(Mbeam+Mtarg); // Velocity of the center of mass
		double Ecm = Mtarg*Beam_E/(Mbeam+Mtarg); // Energy of the center of mass
		double VejectCoM = std::sqrt((2.0/(Meject+Mrecoil))*(Mrecoil/Meject)*(Ecm+Qvalue-(0.0+RecoilExStates[i]))); // Ejectile CoM velocity after reaction
	
		EjectTheta[i] = std::atan2(std::sin(beam_theta_com),(std::cos(beam_theta_com)+(Vcm/VejectCoM))); // Ejectile angle in the lab
		double temp_value = std::sqrt(VejectCoM*VejectCoM-pow(Vcm*std::sin(EjectTheta[i]),2.0));
		double Ejectile_V = Vcm*std::cos(EjectTheta[i]) + temp_value; // Ejectile velocity in the lab frame
		double Ejectile_V2 = Vcm*std::cos(EjectTheta[i]) - temp_value; // Ejectile velocity for double valued solutions
		
		Ejectile_E[i] = 0.5*Meject*Ejectile_V*Ejectile_V; // Ejectile energy in the lab frame
		Recoil_E[i] = (Beam_E+Qvalue-(0.0+RecoilExStates[i])) - Ejectile_E[i]; // Recoil calculations (now in the lab frame)
		RecoilTheta[i] = std::asin(((std::sqrt(2*Meject*Ejectile_E[i]))/(std::sqrt(2*Mrecoil*Recoil_E[i])))*std::sin(EjectTheta[i])); // Recoil angle in the lab
	
		if(VejectCoM >= Vcm){ 
			// Veject is single valued	
			Ejectile_E2[i] = -1;
			Recoil_E2[i] = -1;
			RecoilTheta2[i] = -1;
		} 
		else{ 
			// Veject is double valued
			Ejectile_E2[i] = 0.5*Meject*Ejectile_V2*Ejectile_V2; // Ejectile energy in the lab frame
			Recoil_E2[i] = (Beam_E+Qvalue-(0.0+RecoilExStates[i])) - Ejectile_E2[i]; // Recoil calculations (now in the lab frame)
			RecoilTheta2[i] = std::asin(((std::sqrt(2*Meject*Ejectile_E2[i]))/(std::sqrt(2*Mrecoil*Recoil_E2[i])))*std::sin(EjectTheta[i])); // Recoil angle in the lab		
		}
	

	} // over NrecoilStates
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

// Sample the angular distributions and load them into an array.
// This function assumes values is at least large enough to
// accept a number of doubles equal to the number of distributions
void Kindeux::Sample(double *values){
	if(ang_dist){
		for(unsigned short i = 0; i < NDist; i++){
			values[i] = distributions[i].Sample();
		}
	}
}

// Print debug information
// Does nothing if angular distributions are not set
void Kindeux::Print(){
	if(ang_dist){
		for(unsigned short i = 0; i < NDist; i++){
			std::cout << "  State " << i+1 << ":";
			std::cout << " Reaction X-Section: " << distributions[i].GetReactionXsection() << " mb";
			std::cout << "\tExpected Rate: " << distributions[i].GetRate() << " pps\n";
		}
	}
}
