// rndgauss.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:35:52 2014

double rndgauss0(double w){
	// rndgauss0, a cut down version of rndgauss1; 
	// returns a random number with FWHM w centred at 0;
	
	double w, t, tsq;  
	const double c0=2.515517, c1=0.802853, c2=0.010328; 
	const double d1=1.432788, d2=0.189269, d3=0.001308; 
	const double widthfact=0.424628450; 

	if(w == 0.0){ return 0.0; } 

	t = rand(0); 

	if(t > 0.5){ t = t-0.5; }
	if(t < 1e-30){ t = 11.46380587; }
	else{ 
		tsq = -log(t*t); 
		t = std::sqrt(tsq); 
		
		//     compute inverse by equn 26.2.23
		t=t-(c0+c1*t+c2*tsq)/(1.00+d1*t+(d2+d3*t)*tsq); 
	} 

	//     now randomize x in positive and negative direction
	if(t > 0.5){ t = -t; }
	
	//     now correct for standard deviation
	return widthfact*w*t; 
} 

// rndgauss1 : Generate random numbers with a Gaussian distribution *
void rndgauss1(u,x,f,c,s){ 
	//     a subroutine to calculate one random deviate
	//     for a normal distribution with centroid c and
	//     standard deviation  s
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
	
	double u, x, f, c, s; 

	const double c0=2.5155170, c1=0.8028530, c2=0.0103280); 
	const double d1=1.4327880, d2=0.1892690, d3=0.0013080); 
	const double rcpsqr2pi=0.398942280, sqrt2=1.4142135620); 

	double t, tsq; 
	bool n; 

	u = rand(0); 
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
	//     x=x* (2* nint(drand(0)) -1)
	if(u > 0.5){ x=-x; }

	//     compute function
	f = rcpsqr2pi*std::exp(-(x*x)); 

	//     now correct for centroid and standard deviation
	x = sqrt2*s*x+c; 
	f = f/s; 
} 
