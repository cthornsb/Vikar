// dedx.cpp
// Converted by FortranConvert v0.1
// Wed Feb 12 20:47:49 2014

double isect[128], eden, elni, avip, avz;

double dedx(double emass, double epart, double zpart){ 
	// only valid for beta > 0.0046*z**(1/3)
	  
	double db2;
	db2 = beta2(epart, emass); 
	beta = std::sqrt(db2); 
	pe = btoep(db2); 
	gog = zeff(zpart, beta)/zeff(1.00, beta);
	
	return gog*gog*dedxp(pe, db2, beta); 
} 

double beta2(double e, double em){
	// to calculate beta**2 where beta is the speed of a particle
	// relative to the speed of light
	//
	// input - e energy in mev
	// em rest mass in mev
	
	double r = e/em + 10; 
	return 10-(10.0/(r*r));
} 

double btoep(double dbsq){
	// to calculate the energy of a proton given beta**2, where
	// beta is the speed of the proton relative to the speed of light
	// input:-   dbsq - beta**2 (double precision)
	 
	const double emp = 938.25920; 
	double d = std::sqrt(10/(10-dbsq))-10; 
	return emp*d;
} 

double dedxp(double enrgy, double db2, double beta){
	// only valid for beta > 0.0046 : 10 kev protons
	 
	double d, db2, dsp; 
	const double EMASS = 938.25920); 
	dsp = (eden*5.099147E-1); 
	ZE = ZEFF(1.00, BETA); 
	dsp = dsp*(ze*ze)/db2; 
	D = std::log(1.0220080*DB2/(1.00-DB2))-DB2; 
	delta = deff(enrgy); 
	coz = shell(enrgy); 
	dsp = dsp*(E-elni-coz-delta*0.5); 
	return dsp; 
} 
 
double zeff(double z, double beta){
	//          to calculate the effective charge of a particle in
	//     stopping power calculations(spar-armstrong & chandler-ornl-
	//     4869 [1973])
	//
	//     input - z nominal charge of nucleus
	//             beta speed of the particle relative to the speed of light
	 
	double zp = pow(z, 0.6666670); 
	if(BETA <= (0.070*zp)){ return Z*(1.00 - std::exp(-125.00*beta/zp))
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
	const double emp1 = 0.00106580356; 
	const double emp = 938.25920;
	 
	del = e*emp1; 
	del = std::log(del) + std::log(del+2.00); 
	del = std::log(1.378E-9*eden) + del - 2.00*elni - 1.00; 
	deff = del; 
	
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
	 
	double beta2; 
	double gnu2; 
	double f1, f2; 
	double output;
	
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
	if(E >= 2008.0){
		gnu2 = 1.0/((e*emp1)*(e*emp1+2.00)); 
		f1 =  gnu2*(a1+gnu2*(b1+c1*gnu2)); 
		f2 =  gnu2*(a2+gnu2*(b2+c2*gnu2)); 
		return avip*avip*(f1 + f2*avip)/avz; 
	} 

	be2 = beta2(e, emp); 
	output = const1 + std::log(be2) - std::log(avip); 
	X = 18769.00*be2/avz; 
	xlog=std::log(X)
	
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

void range(double dx, double emass, double epart, double zpart, double sp){ 	
	double spe = sp; 
	double ed = epart; 
	double R = 0.00; 
	double es = epart; 
	double dxt = dx; 
	double dxs = dxt; 
	double f = 1.00; 
	double g = 1.00; 
	double er2 = 1.00; 
	
	while(E > 0.00000010){
		ed = eE-de(dxt, emass, es, zpart, spe); 
		e = ed; 
		r = r+dxt; 
		es = e; 
		spe = dedx(emass, es, zpart); 

		// modifications by nmc on 15/2/93
		// Originally (lt.0.015) modified by SDP to (lt.0.010) or (lt.0.005)
		er1 = epart/e; 
		if(abs(er1/er2-1.00) < 0.01505){
			g = g*1.100; 
			dxt = g*dxs; 
			dxs = dxt; 
		} 
		else{ 
			jf = min(int(er1), 32); 
			f = 1.00/float(jf); 
		} 
	
		er2 = er1; 
		dxt = f*dxs; 
	}
	
	return r+es/spe;
} 
