// test.cpp
// Converted by FortranConvert v0.1
// Mon Feb 17 13:36:09 2014

#include <iostream>
#include <fstream>
#include "vikar_core.h"

int main(){ 
	double r,t,p,x,y,z,b2,btoe; 
	double zf,df,sh,ddxp,ran; 
	double alg,p1,p2,p3,d,lin,mom,rl; 
	double rng,g1,g2,g3,g4,g5,vel,x1,theta,phi; 
	double theta2,phi2,par,v1,v2,v3,v4,ddx; 
	
	ncdedx(1.063,7,3.5,7,4,26.0,p1,p2,p3); 
	std::cout << " " << "ncdedx1:" << " " << p1 << " " << p2 << " " << p3 << " " << std::endl;

	ncdedx(1.063,28,14,7,4,0.5,p1,p2,p3); 
	std::cout << " " << "ncdedx2:" << " " << p1 << " " << p2 << " " << p3 << " " << std::endl;

	ncdedx(1.063,28,14,1,1,0.5,p1,p2,p3); 
	std::cout << " " << "ncdedx3:" << " " << p1 << " " << p2 << " " << p3 << " " << std::endl;

	ran = range(0.001,1,1,1,0.1);
	std::cout << " " << "range1:" << " " << ran << " " << std::endl;

	ran = range(0.02,10,2,3,0.21);
	std::cout << " " << "range2:" << " " << ran << " " << std::endl;
	
	ran = range(0.1,200,4,2,0.01);
	std::cout << " " << "range3:" << " " << ran << " " << std::endl;

	std::cout << " dedx: " << dedx(9.0,12.0,6.0) << std::endl;
	
	std::cout << " dedxp: " << dedxp(26.0,0.25,0.5) << std::endl;
	
	std::cout << " de: " << de(0.01,9.0,12.0,6.0,0.1) << std::endl;
	
	std::cout << " shell: " << shell(26.0) << std::endl;
	
	std::cout << " beta2: " << beta2(2.0,1.0) << std::endl;
	
	std::cout << " zeff: " << zeff(12.0,0.5) << std::endl;
	
	std::cout << " btoep: " << btoep(0.125) << std::endl;
	
	std::cout << " deff: " << deff(10000.0) << std::endl;

	/*std::ofstream temp("new_test.dat");
	for(unsigned short i = 0; i < 10000; i++){
		temp << rndgauss0(4.130) << std::endl;
	}
	temp.close();*/
	
	/*std::ofstream temp("new_test.dat");
	for(unsigned short i = 0; i < 10000; i++){
		strag_targ(12,6,1.320,4.920,1.20,5.60,theta2,phi2,4.920); 
		temp << theta2 << "\t" << phi2 << std::endl;
	}
	temp.close();*/
	
	/*double thetaRecoil, phiRecoil, thetaEject;
	double phiEject, Erecoil, Eeject;
	
	Kindeux kind;
	kind.Initialize(7.0,2.0,8.0,1.0,-2.0);
	std::ofstream temp("new_test.dat");
	for(unsigned short i = 0; i < 10000; i++){
		kind.FillVars(26.0,0.0,3.14,0.0,0.0,thetaRecoil,phiRecoil,thetaEject,phiEject,Erecoil,Eeject);
		temp << phiEject << "\t" << Eeject << "\n";
	}
	temp.close();*/
	
	return 0;
} 
