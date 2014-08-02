// vikar_lib.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 19:33:21 2014

#include "vikar_core.h"
#include "planar.h"

/////////////////////////////////////////////////////////////////////
// Constant Globals (for fortran commons)
/////////////////////////////////////////////////////////////////////

const double pi = 3.1415926540;
const double deg2rad = pi/180.0;
const double rad2deg = 180.0/pi;
double eden = 0.0;
double elni = 0.0;
double avip = 0.0;
double avz = 0.0;
double E = 0.0;

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

/////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////

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
template <typename T>
std::string to_str(T input){
	std::stringstream output;
	output << input;
	return output.str();
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
	//std::cout << tgtdens << "," << atarget << "," << ztarget << "," << abeam << "," << zbeam << "," << energy << std::endl;	
	double amun[4], xsn[4];
	double mchem[10], zmed[10], amua[10];
	double nmed = 1.0;
	double avden = tgtdens;
	amua[0] = atarget;
	zmed[0] = ztarget;
	double press = 760.0;
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
	
	double t, tsq;  
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
	
	// Dignostics - to be removed whan satisfied
	//      energy = 5.80
	//      Z = 1.0 !4
	//      A = 1.0
	//      thickness = 0.2 !25.000
	//
	//      v = sqrt(2.0*energy/(A*931.5))
	//      p = sqrt(2.0*A*931.5*energy)
	
	// Most of this (the v and p calcualtion) cancels out in the equation
	// It's left in for clarity, but I might remove it, as it's pointless
	double v = velocity(energy,A); 
	double p = momentum(energy,A); 
	
	// Dignostic - to be removed whan satisfied
	theta = 13.6/(v*p)*Z*std::sqrt(thickness/X)*(1.0+0.038*std::log(thickness/X)); 
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
	
	double norm; 
	
	// Put in check of length > 0    ?????
	length = std::sqrt(pow(xl, 2)+pow(yl, 2)+pow(zl, 2)); 
	norm = std::sqrt(1.0/(pow(xl, 2)+pow(yl, 2)+pow(zl, 2))); 
	x = xl*norm; 
	y = yl*norm; 
	z = zl*norm; 
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
    	
void Kindeux::Initialize(double Mbeam_, double Mtarg_, double Mrecoil_, double Meject_, double Qvalue_){
	Mbeam = Mbeam_;
	Mtarg = Mtarg_;
	Mrecoil = Mrecoil_;
	Meject = Meject_;
	Qvalue = Qvalue_;
}

void Kindeux::SetDist(unsigned short NDistPoints_, double Int_max_, double *DistAng_, double *DistInt_){
	NDistPoints = NDistPoints_;
	Int_max = Int_max_;
	DistAng = DistAng_;
	DistInt = DistInt_;
	ang_dist = true;
}

void Kindeux::FillVars(double Ebeam, double theta_beam, double phi_beam, double ExEject, double ExRecoil, double &RecoilTheta, 
	      double &RecoilPhi, double &EjectTheta, double &EjectPhi, double &Erecoil, double &Eeject){
	// kindeux 1.0 written by S.D.Pain on 24/11/2004
	//
	// kindeux 2.0 updated by S.D.Pain on 4/3/2005
	//    - Updated to utilise transform3
	//    - Cleaned up old junk, diagnostics etc
	//
	// kindeux 3.0 updated by S.D.Pain on 5/5/2006
	//    - Updated to employ non-isotropic angular distributions
	//
	// Subroutine for calculating two body reaction kinematics
	// The reaction energy is read in, along with the direction of the incident
	// ion. The Q value, excitation energies and particle masses are read in, also.
	// The reaction occurs evenly distributed over all solid angle, in the CoM frame,
	// or is defined by an angular distribution profile for the ejectile, passed
	// from the parent program. This profile is in the form of a cumulative integral
	// of an angular distribution of the form dSigma/dTheta (not dSigma/dOmega) in
	// the CoM system, with 0 degrees defined by normal kinematics.
	// The lab energies and directions are calculated for the recoil and ejectile.
	// The recoil and ejectile directions are altered due to the angle of the
	// incoming beam particle. E, theta, phi are returned for recoil and ejectile.

	// Internal Variables
	double E, vBeam, vCoM, vRecoil, vEject; 
	double recoilX, recoilY, recoilZ; 
	double ejectX, ejectY, ejectZ; 
	double beamX, beamY, beamZ; // Lab direction of beam
	double CoMX, CoMY, CoMZ; // CoM direction of beam

	// Calculate beam velocity
	vBeam = velocity(Ebeam,Mbeam);
	 
	// Calculate the CoM velocity in the lab
	vCoM = vBeam*Mbeam/(Mbeam+Mtarg); 

	// Set the direction of the CoM frame to along the Z axis.
	// Perform all calcualtions in this frame, and rotate the
	// lab vectors later.

	// Convert the CoM velocity to cartesions (should be entirely along
	// the Z direction), and convert to a unit vector. Proabably
	// pointless now, as everything is aligned along the Z axis
	Sphere2Cart(vCoM,0.0,0.0,beamX,beamY,beamZ); 
	unitV(beamX,beamY,beamZ,CoMX,CoMY,CoMZ,vCoM); 

	// Randomly select outgoing coordinates in CoM frame
	// I think this samples evenly over solid angle... check
	// for the ejectile
	// NB theta = 0 corresponds to forward angles in the lab
	// so for the inverse kinematics case, need to reverse the
	// coordinates so theta = 0 gives backward angles in the lab
	EjectPhi = 2.0*pi*frand(); 

	if(ang_dist){
		//i = 1; 
		//while(temp > DistInt[i]){ i = i+1; } 
		double temp;
		unsigned short i;
		temp = frand()*Int_max;
		for(i = 0; i < NDistPoints; i++){
			if(temp <= DistInt[i]){ break; }
		}
	
		EjectTheta = linear(DistInt[i-1],DistInt[i],DistAng[i-1],DistAng[i],temp); 
		EjectTheta = EjectTheta*deg2rad; 
	
		// swap directions for inverse kinematics cases
		if(Mbeam > Mtarg){ EjectTheta = pi-EjectTheta; }
	} 
	else{ EjectTheta = acos(-2.0*frand()+1.0); } 

	// Calcualte the recoil's theta and phi
	RecoilTheta = pi-EjectTheta; 
	if(EjectPhi < pi){ RecoilPhi = EjectPhi+pi; }
	else{ RecoilPhi = EjectPhi-pi; }
	//if(EjectPhi >= pi){ RecoilPhi = EjectPhi-pi; }

	// Calculate the CoM energy, and add on the Q value and subtract
	// any excitations
	E = 0.5*Mbeam*pow((vBeam-vCoM), 2.0) + 0.5*Mtarg*pow(vCoM, 2.0); 
	E += (Qvalue - ExEject - ExRecoil);

	// Calculate the magnitudes of the velocities of the reaction products
	// in the CoM frame
	vRecoil = std::sqrt((2.0*E)/(Mrecoil+pow(Mrecoil, 2.0)/Meject)); 
	vEject = std::sqrt((2.0*E)/(Meject+pow(Meject, 2.0)/Mrecoil)); 

	// Convert the recoil's and ejectile's CoM velocity vectors from
	// spherical polars to cartesians, and transform into the laboratory frame
	Sphere2Cart(vRecoil,RecoilTheta,RecoilPhi,recoilX,recoilY,recoilZ); 
	Sphere2Cart(vEject,EjectTheta,EjectPhi,ejectX,ejectY,ejectZ); 

	recoilX += vCoM*CoMX; 
	recoilY += vCoM*CoMY; 
	recoilZ += vCoM*CoMZ; 
	ejectX += vCoM*CoMX; 
	ejectY += vCoM*CoMY; 
	ejectZ += vCoM*CoMZ; 

	// Convert the lab frame cartesians to lab frame spherical polars
	// NB now theta_* and phi_* contain lab angles, not CoM angles
	// These variables are returned to the Master
	Cart2Sphere(recoilX,recoilY,recoilZ,vRecoil,RecoilTheta,RecoilPhi); 
	Cart2Sphere(ejectX,ejectY,ejectZ,vEject,EjectTheta,EjectPhi);

	// Calculate the laboratory energy
	Erecoil = 0.5*Mrecoil*pow(vRecoil, 2.0); 
	Eeject = 0.5*Meject*pow(vEject, 2.0); 

	// Rotate the velocity vectors due to the incident angle of the beam particle
	transform(theta_beam,phi_beam,EjectTheta,EjectPhi,EjectTheta,EjectPhi); 	
	transform(theta_beam,phi_beam,RecoilTheta,RecoilPhi,RecoilTheta,RecoilPhi);
}

/////////////////////////////////////////////////////////////////////
// product_proc.f
/////////////////////////////////////////////////////////////////////

void ProdProc::Initialize(double A_, double Z_, double beamspot_, double thickness_, double targ_depth_, double targ_angle_, double *E_targ_, double *Erange_targ_, 
			  double targRadL_, unsigned short Ndet_plan_, unsigned short *Nstrips_plan_, double *DetZ_plan_, double **DetXMax_, double **DetXMin_, 
			  double **DetYMax_, double **DetYMin_, double *E_det_, double *Erange_det_, double conv_det_, double DetRadL_, double *dEthick_plan_, 
			  double *Ethick_plan_, double *Eres_dE_plan_, double *Pres_dE_plan_, double *Eres_E_plan_, double *Pres_E_plan_, double *Pres_dE_en_plan_, 
			  double *Pres_E_en_plan_, bool *xPlane_, bool *yPlane_, double *detWidth_plan_, double *detLength_plan_){	
	A = A_; 
	Z = Z_;
	beamspot = beamspot_;
	thickness = thickness_;
	targ_depth = targ_depth_;
	targ_angle = targ_angle_;
	targRadL = targRadL_;
	E_targ = E_targ_;
	Erange_targ = Erange_targ_;
	Ndet_plan = Ndet_plan_;
	Nstrips_plan = Nstrips_plan_;
	DetZ_plan = DetZ_plan_;
	E_det = E_det_;
	Erange_det = Erange_det_;
	dEthick_plan = dEthick_plan_;
	Ethick_plan = Ethick_plan_;
	DetXMax = DetXMax_;
	DetXMin = DetXMin_;
	DetYMax = DetYMax_;
	DetYMin = DetYMin_;
	conv_det = conv_det_;
	DetRadL = DetRadL_;
	Eres_dE_plan = Eres_dE_plan_;
	Pres_dE_plan = Pres_dE_plan_;
	Eres_E_plan = Eres_E_plan_;
	Pres_E_plan = Pres_E_plan_;
	Pres_dE_en_plan = Pres_dE_en_plan_;
	Pres_E_en_plan = Pres_E_en_plan_;
	detWidth_plan = detWidth_plan_;
	detLength_plan = detLength_plan_;
	xPlane = xPlane_;
	yPlane = yPlane_;
	init = true;	
}

// Return true if hit is detected and false otherwise
bool ProdProc::FillVars(double energy, double theta, double phi, unsigned short &DetHitPlan, unsigned short &StripHitPlan, double &E_emerge,
		  double &theta_Emerge, double &phi_Emerge, double &det_dE, double &det_E, double &theta_E, double &phi_E, double &det_dE_final,
		  double &det_E_final, double &theta_dE_final, double &phi_dE_final, double &theta_E_final, double &phi_E_final, bool &stuck){
	if(!init){ return false; }
	// product_proc 1.0 written by S.D.Pain in the depths of history
	//
	// product_proc 1.1 updated by S.D.Pain on 11/12/2013
	// to improve handling of beamspot effects for annular detectors (effect previously over-estimated)
	//
	// product_proc 1.1 VANDLE updated by S.D.Pain on 2014/01/24
	// to remove energy loss effects for neutrons (three lines, marked SDP, VANLDE, neutrons etc
	//
	// This subroutine processes the life of a reaction product after
	// it is emitted from the reaction. This includes the effect of
	// passing out of the target, and all detector effects
	//
	// Note that a detector threshold of 0.2 MeV has been applied
	// Note: detector threshold changed to 0.1 MeV by Cory (for small bars)
	
	// Basic variables
	double theta_Emerge_old; 
	
	// Target variables
	double range_targ, depth_out; 
	double theta_new, phi_new; 
	double lost_E;
	
	// Detector variables
	unsigned short detNo, stripNo; 
	
	// Detector Energy loss variables
	double range_det; 
	double DetThick; 
	
	// Detector Resolution Variables
	unsigned short i; 
	
	// Calculate the particle's range in the target
	i = 0; 
	while(E_targ[i] < energy){ i = i + 1; } 
	range_targ = linear(E_targ[i-1],E_targ[i],Erange_targ[i-1],Erange_targ[i],energy); 

	// Caculate the target thickness (depth_out) the particle sees
	targ_thick(theta,phi,thickness,targ_depth,targ_angle,depth_out); 

	// If the particle stops in the target, exit the event
	if(depth_out >= range_targ){ // !!!!!!!!!!!!!!!! Fix this!!!
		stuck = true;
		return false;
	} 

	// Find the ranges for energies straddling the particle's energy
	i = 0; 
	while(Erange_targ[i]<(range_targ-depth_out)){
		i = i + 1; 
	} 
	
	// calculate the energy of the emergant particle
	E_emerge = linear(Erange_targ[i-1],Erange_targ[i],E_targ[i-1],E_targ[i],(range_targ-depth_out)); 
	E_emerge = energy; // SDP - for neutrons,VANDLE,Bill!

	// Angular Straggling of ejectile in the target
	strag_targ(A,Z,depth_out,theta,phi,energy,theta_new,phi_new,targRadL);
	theta_Emerge = theta_new; 
	phi_Emerge = phi_new; 

	theta_Emerge = theta; // SDP - for neutrons,VANDLE,Bill!
	phi_Emerge = phi; // SDP - for neutrons,VANDLE,Bill!

	// End of passage through target processing
	//---------------------------------------------------------------------------

	//--------------------------------------------------------------------------- 
	// Detection system processing
	unsigned short HitCounter = 0; 
	if(E_emerge > 0.2 && E_emerge < 20.0){ // SDP - normally 0.5 - Typ. 0.1 to 0.2 MeV for v-bars - Cory
		// See if it hits a planar detector
		if(Ndet_plan > 0){
			DetHit_plan(DetXMax,DetXMin,DetYMax,DetYMin,Ndet_plan,Nstrips_plan,DetZ_plan,theta_Emerge,phi_Emerge,HitCounter,DetHitPlan,StripHitPlan); 
		} 

		theta_Emerge_old = theta_Emerge; // SDP beamspot size adjustment, 2013-10-15
	} // threshold check

	theta_Emerge = theta_Emerge_old; // SDP beamspot size adjustment, 2013-10-15
	
	// A Good hit has been found. Continue with detection
	if(HitCounter == 1){
		//------------------------------------------------------------------; 
		// dE DETECTOR
		// Calculate the ejectile's range in the dE detector material
		i = 0; 
		while(E_det[i] < E_emerge){ i = i + 1; } 
		range_det = linear(E_det[i-1],E_det[i],Erange_det[i-1],Erange_det[i],E_emerge); 

		// Calculate the thickness the ejectile sees in the dE layer
		det_thick_plan(theta_Emerge,phi_Emerge,dEthick_plan[DetHitPlan],DetThick); 

		// If the particle stops in the dE detector
		if(DetThick >= range_det){
			det_dE = E_emerge; 
			det_E = 0.0; 
		} 
		else{ 
			// Find the ranges for energies straddling the ejectile energy		 
			i = 0; 
			while(Erange_det[i] < (range_det-DetThick)){
				i = i + 1; 
			} 

			// calculate the energy of the emergant ejectile
			det_E = linear(Erange_det[i-1],Erange_det[i],E_det[i-1],E_det[i],(range_det-DetThick)); 
			det_dE = E_emerge - det_E; 

			// **************************************************
			// Angular Straggling
			 strag_dE_plan(A,Z,DetThick,DetZ_plan[DetHitPlan],theta_Emerge,phi_Emerge,E_emerge,theta_E,phi_E,DetRadL,conv_det); 
		} // Particle stops in dE check

		//------------------------------------------------------------------; 
		// E DETECTOR
		// Calculate the ejectile's range in the E detector material
		if(det_E > 0){// SDP 1.5 MeV threshold
			i = 0; 
			while(E_det[i] < det_E){ i = i + 1; } 
			range_det = linear(E_det[i-1],E_det[i],Erange_det[i-1],Erange_det[i],det_E); 

			// Calculate the thickness the ejectile sees in the E layer
			//  Note that this uses the effective straggling angle rather than the real angle!!!
			//  Correct this by getting AngStrag to return the real angle as well!!
			det_thick_plan(theta_E,phi_E,Ethick_plan[DetHitPlan],DetThick); 

			// If the particle punches through the E detector
			if(DetThick<range_det){
				// Find the ranges for energies straddling the ejectile energy
				i = 0; 
				while(Erange_det[i] < (range_det-DetThick)){ i = i + 1; } 
			
				// calculate the energy of the emergant ejectile
				lost_E = linear(Erange_det[i-1],Erange_det[i],E_det[i-1],E_det[i],(range_det-DetThick)); 
				det_E = det_E - lost_E; 
			} // Particle stops in E check
		} // Particle gets to E check

		// End of E DETECTOR
		//------------------------------------------------------------------; 

		//------------------------------------------------------------------; 
		// RESOLUTIONS
		// Apply the detector resolutions to the particles, and calcualte the
		// apparent energies and angles
		det_dE_final = det_dE; 
		theta_dE_final = theta_Emerge; 
		det_E_final = det_E; 
		theta_E_final = theta_Emerge; 

		resolution_plan(det_dE,theta_Emerge,phi_Emerge,Eres_dE_plan[DetHitPlan],Pres_dE_plan[DetHitPlan],Pres_dE_en_plan[DetHitPlan],
				DetZ_plan[DetHitPlan],xPlane[DetHitPlan],yPlane[DetHitPlan],DetXMin[DetHitPlan][StripHitPlan],DetYMin[DetHitPlan][StripHitPlan],
				detWidth_plan[DetHitPlan],detLength_plan[DetHitPlan],Nstrips_plan[DetHitPlan],StripHitPlan,det_dE_final,theta_dE_final,phi_dE_final); 
			
		if(det_E == 0.0){
			det_E_final = 0.0; 
			theta_E_final = 0.0; 
		} 
		else{ 
			resolution_plan(det_E,theta_E,phi_E,Eres_E_plan[DetHitPlan],Pres_E_plan[DetHitPlan],Pres_E_en_plan[DetHitPlan],
					DetZ_plan[DetHitPlan],xPlane[DetHitPlan],yPlane[DetHitPlan],DetXMin[DetHitPlan][StripHitPlan],DetYMin[DetHitPlan][StripHitPlan],
					detWidth_plan[DetHitPlan],detLength_plan[DetHitPlan],Nstrips_plan[DetHitPlan],StripHitPlan,det_E_final,theta_E_final,phi_E_final); 
		} 
		
		return true;
	} // HitCounter = 1 
	else if(HitCounter > 1){ 
		std::cout << "blair!\n"; 
		return false;
	}
	
	return false;
};
