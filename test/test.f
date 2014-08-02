      program vikar
      implicit none
      
        double precision r,t,p,x,y,z,b2,btoe
      	double precision zf,df,sh,ddxp,ran
        double precision alg,p1,p2,p3,d,lin,mom,rl
        double precision rng,g1,g2,g3,g4,g5,vel,x1,theta,phi
        double precision x2,y2,par,v1,v2,v3,v4
        integer i
            	
        double precision beta2,btoep,zeff,deff,shell,dedxp
        double precision range2,algip,de,linear,momentum
        double precision radlength,rndgauss0,velocity,ddx,dedx
      	double precision thetaRecoil,phiRecoil,thetaEject
      	double precision phiEject,Erecoil,Eeject
      
      	double precision NAngPoints,maxInt,DistAng_temp
      	double precision DistInt_temp,ADist_type

      	call ncdedx(1.063d0,7d0,3.5d0,7d0,4d0,26.0d0,p1,p2,p3)
      	write(6,*)" ncdedx:",p1,p2,p3

      	call ncdedx(1.063d0,28d0,14d0,7d0,4d0,0.5d0,p1,p2,p3)
      	write(6,*)" ncdedx2:",p1,p2,p3

	call ncdedx(1.063d0,28d0,14d0,1d0,1d0,0.5d0,p1,p2,p3); 
	write(6,*)" ncdedx3:",p1,p2,p3
      	
      	ran = range2(0.001d0,1d0,1d0,1d0,0.1d0);
      	write(6,*)" range:",ran
      	
	ran = range2(0.02d0,10d0,2d0,3d0,0.21d0);
	write(6,*)" range2:",ran
	
	ran = range2(0.1d0,200d0,4d0,2d0,0.01d0);
	write(6,*)" range3:",ran
      	
      	write(6,*)" dedx: ",dedx(9.0d0,12.0d0,6.0d0)
      	
      	write(6,*)" dedxp: ",dedxp(26.0d0,0.25d0,0.5d0)
      	
      	write(6,*)" de: ",de(0.01d0,9.0d0,12.0d0,6.0d0,0.1d0)
      	
      	write(6,*)" shell: ",shell(26.0d0)
      	
      	write(6,*)" beta2: ",beta2(2.0d0,1.0d0)
      	
      	write(6,*)" zeff: ",zeff(12.0d0,0.5d0)
      	
      	write(6,*)" btoep: ",btoep(0.125d0)
      	
      	write(6,*)" deff: ",deff(10000.0d0)

c      	open(unit=11,file="old_test.dat",status="unknown")
c      	do i=1,10000
c      	      	write(11,*)rndgauss0(4.13d0)
c      	enddo
c      	close(11)

c      	open(unit=11,file="old_test.dat",status="unknown")
c      	do i=1,10000
c      	      	call strag_targ(12d0,6d0,1.32d0,4.92d0,1.2d0,5.6d0,x2,y2,4.92d0)
c      	      	write(11,*)x2,y2
c      	enddo

c      	do i=1,10000
c          call kindeux(26d0,0d0,3.14d0,-2d0,0d0,0d0,
c     & 7d0,2d0,8d0,1d0,ADist_type,
c     & NAngPoints,maxInt,DistAng_temp,
c     & DistInt_temp,thetaRecoil,phiRecoil,thetaEject,phiEject,
c     & Erecoil,Eeject)!returned data

c          write(11,*) phiEject,Eeject ! (CORY)
c      	enddo
c      	close(11)

      end
