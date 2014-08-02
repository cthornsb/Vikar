// ncdedx.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:35:04 2014

// DEDX : Subroutine to calculate energy loss of ions in a medium * 
void ncdedx(tgtdens,atarget,ztarget,abeam,zbeam,energy,dedxmg,tgtionpot,rangemg){ 
	common isect(128),eden,elni,avip,avz; 
	dimension amun(4),xsn(4); 
	dimension mchem(10),zmed(10),amua(10); 
	external range; 
	//
	//     mchem(i) - # atoms of element(i) in molecule of medium
	//     zmed(i)  - atomic number of element(i)
	//     amua(i)  - mass of one atom of element(i) in a.m.u.
	//
	nmed=1; 
	AVDEN=tgtdens; 
	AMUA(1)=atarget; 
	ZMED(1)=ztarget; 
	press=760.0; 
	mchem(1)=1; 
	amum=mchem(1)*amua(1); 
	denm=avden/(amum*1.660543); 
	//
	//     calculate electrons/cc, av. ionization potl, av. charge
	//
	SUMN=0.00; 
	SUMNZ=0.00; 
	SUMNZI=0.00; 
	do 120 i=1,nmed; 
	en=mchem(i)*denm; 
	enz=en*zmed(i); 
	enzi=enz*algip(zmed(i)); 
	sumn=sumn+en; 
	sumnz=sumnz+enz; 
	120 sumnzi=sumnzi+enzi; 
	elni=sumnzi/sumnz; 
	avip=exp(elni); 
	tgtionpot=AVIP*1000000; 
	eden=sumnz; 
	avz=sumnz/sumn; 
	//
	//          read in parameters of beam to be analyzed
	//
	nnuc=1; 
	AMUN(1)=abeam; 
	ZPART=zbeam; 
	102 xsn(1)=0.0; 
	210 AMASS=AMUN(1)+XSN(1)/931481.20; 
	EMASS=AMASS*931.48120; 
	//
	//          print out details of medium
	//
	EPART = energy; 
	// Set line below from 0.5/() to 0.05/() to remove discontinuities
	// in calculated ranges. Not checked validiy yet.
	// SDP 16/12/2004
	DX = 0.050/DEDX(EMASS,EPART,ZPART); 
	226 sp=dedx(emass,epart,zpart); 
	r=range(dx,emass,epart,zpart,sp); 
	
	dedxmg=(SP*0.001)/AVDEN; 
	rangemg=R*AVDEN*1000.00; 
	
	return 0; 
} 
 
double algip(double z){ 
	//          to calculate alog(ionization potential) for an element
	//     of atomic number z
	//
	//      n.b. ionization potl in mev !!!!!!
	iz = z + 0.050; 
	if(iz > 13){ potl = 9.760*z + 58.80/(pow(z, 0.190)); }
	else{ potl=pot(iz); }
	
	return std::log(potl * 1.0e-6); 
} 
