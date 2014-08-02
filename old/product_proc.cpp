// product_proc.cpp
// Converted by FortranConvert v0.1
// Thu Feb 13 13:35:26 2014

#ifndef PRODUCT_PROC_CPP
#define PRODUCT_PROC_CPP

#include <cmath>
#include <iostream>

#include "strag_dE_plan.cpp"
#include "Det_Planar.cpp"
#include "Res_Planar.cpp"

void product_proc(double energy, double theta, double phi, double A, double Z, double beamspot, double thickness, double targ_depth, double targ_angle, double *E_targ, double *Erange_targ, double targRadL,
		  int NDet_cyl, int *Nstrips_cyl, double **DetPhiMin_cyl, double **DetPhiMax_cyl, double **DetThetaMin_cyl, double **DetThetaMax_cyl, double *detLength_cyl, double *detWidth_cyl, double *DetRad_cyl, double *DetZ_cyl, double *DetPhi_cyl, 
		  int Ndet_plan, int *Nstrips_plan, double *DetZ_plan, double **DetXMax, double **DetXMin, double **DetYMax, double **DetYMin, 
		  int Ndet_ann, int *Nstrips_ann, double **DetPhiMin_ann, double **DetPhiMax_ann, double **DetThetaMin_ann, double **DetThetaMax_ann,
		  double *E_det, double *Erange_det, double conv_det, double DetRadL, 
		  double *dEthick_cyl, double *dEthick_plan, double *dEthick_ann, double *Ethick_cyl,  double *Ethick_plan, double *Ethick_ann,
		  double *Eres_dE_cyl, double *Pres_dE_cyl, double *Eres_E_cyl, double *Pres_E_cyl, double *Pres_dE_en_cyl, double *Pres_E_en_cyl,
		  double *Eres_dE_plan, double *Pres_dE_plan, double *Eres_E_plan, double *Pres_E_plan, double *Pres_dE_en_plan, double *Pres_E_en_plan, 
		  bool *xPlane, bool *yPlane, double *detWidth_plan, double *detLength_plan, double *Eres_dE_ann, double *Eres_E_ann, double *DetZ_Ann,
		  double *DetInner_ann, double *DetOuter_ann, double *DetPhi_ann,
		  int &HitCounter, int &HitGeo, int &DetHitCyl, int &StripHitCyl, int &DetHitPlan, int &StripHitPlan, int &DetHitAnn, int &StripHitAnn,
		  double &E_emerge, double &theta_Emerge, double &phi_Emerge, double &det_dE, double &det_E, double &theta_E, double &phi_E, double &det_dE_final, double &det_E_final, 
		  double &theta_dE_final, double &phi_dE_final, double &theta_E_final, double &phi_E_final, bool &stuck){
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
	double d_radius, d_theta, theta_Emerge_old; 
	
	// Target variables
	double range_targ, depth_out; 
	double theta_new, phi_new; 
	double rad_temp, ran_theta, lost_E;
	
	// Detector variables
	int detNo, stripNo; 
	
	// Detector Energy loss variables
	double range_det; 
	double DetThick; 
	
	// Detector Resolution Variables
	int i; 
	
	HitGeo = 0; 
	stuck = false; 
	
	//---------------------------------------------------------------------------; 
	// Passage through target processing
	//
	//
	//
	//        std::cout << 'here too!'
	//          call resolution_plan(energy,theta,Phi,
	//     &     Eres_dE_plan[DetHitPlan],
	//     &     Pres_dE_plan[DetHitPlan],DetZ_plan[DetHitPlan],
	//     &     detWidth_plan[DetHitPlan],
	//     &     xPlane[DetHitPlan],yPlane[DetHitPlan],
	//     &     DetxMin[DetHitPlan][StripHitPlan],
	//     &     DetyMin[DetHitPlan][StripHitPlan],
	//     &     Nstrips_plan[DetHitPlan],StripHitPlan,
	//     &     dE_final_eject,dE_final_Theta,dE_final_Phi)
	
	// Calculate the particle's range in the target
	i = 0; 
	while(E_targ[i]<energy){ i = i + 1; } 
	range_targ = linear(E_targ[i-1],E_targ[i],Erange_targ[i-1],Erange_targ[i],energy); 

	// Caculate the target thickness (depth_out) the particle sees
	targ_thick(theta,phi,thickness,targ_depth,targ_angle,depth_out); 

	// If the particle stops in the target, exit the event
	if(depth_out>=range_targ){ // !!!!!!!!!!!!!!!! Fix this!!!
		stuck = true;
		return;
	} 

	// Find the ranges for energies straddling the particle's energy
	i = 0; 
	while(Erange_targ[i]<(range_targ-depth_out)){
		i = i + 1; 
		//if(i > 100000){ pause; }
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
	HitCounter=0; 
	
	// Work out (roughly) if a detector was hit
	// For cylindrical polar detectors...
	if(E_emerge>0.2){ // SDP - normally 0.5 - Typ. 0.1 to 0.2 MeV for v-bars - Cory
		if(NDet_cyl>0){
			for(detNo = 0; detNo < NDet_cyl; detNo++){
				for(stripNo = 1; stripNo < Nstrips_cyl[detNo]; stripNo++){
					if(DetPhiMin_cyl[detNo][stripNo]>DetPhiMax_cyl[detNo][stripNo]){
						if((phi_Emerge>=DetPhiMin_cyl[detNo][stripNo] &&phi_Emerge<=(2*pi)) || (phi_Emerge>=(0.0) &&phi_Emerge<=DetPhiMax_cyl[detNo][stripNo])){
							if(theta_Emerge>=DetThetaMin_cyl[detNo][stripNo] &&theta_Emerge<=DetThetaMax_cyl[detNo][stripNo]){
								HitCounter = HitCounter + 1; 
								DetHitCyl = detNo; 
								StripHitCyl = stripNo; 
							} 
						} 
					} 
					else{ 
						if(phi_Emerge >= DetPhiMin_cyl[detNo][stripNo] && phi_Emerge <= DetPhiMax_cyl[detNo][stripNo]){
							if(theta_Emerge>=DetThetaMin_cyl[detNo][stripNo] &&theta_Emerge<=DetThetaMax_cyl[detNo][stripNo]){
								HitCounter = HitCounter + 1; 
								DetHitCyl = detNo; 
								StripHitCyl = stripNo; 
							} 
						} 
					} 
				} 
			} 
		} 

		// If it was hit by a rough calculation, check precisely
		HitGeo = 0; 

		// See if it hits a cylindrical detector
		/*if(HitCounter==1){
			DetHit_cyl(detLength_cyl[DetHitCyl],detWidth_cyl[DetHitCyl],DetRad_cyl[DetHitCyl],DetZ_cyl[DetHitCyl],DetPhi_cyl[DetHitCyl],theta_Emerge,phi_Emerge,HitCounter); 
		} 

		if(HitCounter==1){ HitGeo = 1; }*/

		// See if it hits a planar detector
		if(HitCounter==0){
			if(Ndet_plan > 0){
				DetHit_plan(DetXMax,DetXMin,DetYMax,DetYMin,Ndet_plan,Nstrips_plan,DetZ_plan,theta_Emerge,phi_Emerge,HitCounter,DetHitPlan,StripHitPlan); 
				if(HitCounter==1){ HitGeo = 2; }
			} 
		} 

		theta_Emerge_old = theta_Emerge; // SDP beamspot size adjustment, 2013-10-15

		// Annular hit checks
		/*if(HitCounter==0){
			if(Ndet_ann>0){
				d_radius = rndgauss0(beamspot); // SDP beamspot size adjustment, 2013-12-11
				d_radius = abs(d_radius); // SDP beamspot size adjustment, 2013-12-11
				ran_theta = frand()*pi; // SDP beamspot size adjustment, 2013-12-11
				d_radius = d_radius*std::cos(ran_theta); // SDP beamspot size adjustment, 2013-12-11
				theta_Emerge = theta_Emerge_old; // SDP beamspot size adjustment, 2013-10-15

				for(detNo = 0; detNo < Ndet_ann; detNo++){
					// d_theta = atan(d_radius/DetZ_Ann[detNo]); // SDP beamspot size adjustment, 2013-10-15

					rad_temp = DetZ_Ann[detNo]*tan(theta_Emerge_old); // SDP beamspot size adjustment, 2013-12-11
					rad_temp = rad_temp+d_radius; // SDP beamspot size adjustment, 2013-12-11
					theta_Emerge = atan(rad_temp/DetZ_Ann[detNo]); // SDP beamspot size adjustment, 2013-12-11
					if(DetZ_Ann[detNo]<0){ theta_Emerge = pi+theta_Emerge; } // SDP beamspot size correction for back hemisphere, 2013-12-1

					for(stripNo = 0; stripNo < Nstrips_ann[detNo]; stripNo++){
						if(DetPhiMin_ann[detNo][stripNo]>DetPhiMax_ann[detNo][stripNo]){
							if((phi_Emerge>=DetPhiMin_ann[detNo][stripNo] &&phi_Emerge<=(2*pi)) || (phi_Emerge>=(0.0) &&phi_Emerge<=DetPhiMax_ann[detNo][stripNo])){
								if(theta_Emerge>=DetThetaMin_ann[detNo][stripNo] &&theta_Emerge<=DetThetaMax_ann[detNo][stripNo]){
									HitCounter = HitCounter + 1; 
									DetHitAnn = detNo; 
									StripHitAnn = stripNo; 
								} 
							} 
						} 
						else{ // Det does overlap with Phi=0
							if(phi_Emerge>=DetPhiMin_ann[detNo][stripNo] &&phi_Emerge<=DetPhiMax_ann[detNo][stripNo]){
								if(theta_Emerge>=DetThetaMin_ann[detNo][stripNo] &&theta_Emerge<=DetThetaMax_ann[detNo][stripNo]){
									HitCounter = HitCounter + 1; 
									DetHitAnn = detNo; 
									StripHitAnn = stripNo; 
								} 
							} 
						} 
					} 
				} 
				
				if(HitCounter==1){ HitGeo = 3; }
			} 
		}*/
	} // threshold check

	theta_Emerge = theta_Emerge_old; // SDP beamspot size adjustment, 2013-10-15

	if(HitCounter>1){ 
		std::cout << "blair!\n"; 
		std::cout << NDet_cyl;
		std::cout << "double!\n";
		return;
	}
	
	// A Good hit has been found. Continue with detection
	if(HitCounter==1){
		//------------------------------------------------------------------; 
		// dE DETECTOR
		// Calculate the ejectile's range in the dE detector material
		i = 0; 
		while(E_det[i] < E_emerge){ i = i + 1; } 
		range_det = linear(E_det[i-1],E_det[i],Erange_det[i-1],Erange_det[i],E_emerge); 

		//         if(E_emerge.gt.69.0 .and. E_emerge.lt.91.0 )then
		// std::cout << E_emerge,range_det; 
		// write(58,*) E_emerge,range_det; 
		// endif

		// Calculate the thickness the ejectile sees in the dE layer
		if(HitGeo==1){ 
			//det_thick_cyl(theta_Emerge,phi_Emerge,DetPhi_cyl(detHitCyl),dEthick_cyl(detHitCyl),DetThick);
			std::cout << " Warning: Cylinder hit: " << HitCounter << " " << HitGeo << std::endl; 
		}
		else if(HitGeo==2){ 
			det_thick_plan(theta_Emerge,phi_Emerge,dEthick_plan[DetHitPlan],DetThick); 
		} 
		else if(HitGeo==3){ 
			//det_thick_ann(theta_Emerge,phi_Emerge,dEthick_ann(detHitAnn),DetThick); 
			std::cout << " Warning: Annular hit: " << HitCounter << " " << HitGeo << std::endl; 
		}

		// If the particle stops in the dE detector
		if(DetThick >= range_det){
			det_dE = E_emerge; 
			det_E = 0.0; 
		} 
		else{ 
			// Find the ranges for energies straddling the ejectile energy
			//          if(E_emerge.gt.69.0 .and. E_emerge.lt.91.0 )then
			// std::cout << E_emerge,range_det,DetThick; 
			// write(57,*) E_emerge,range_det,DetThick; 
			// endif
			 
			i = 0; 
			while(Erange_det[i] < (range_det-DetThick)){
				i = i + 1; 
				// if(E_emerge.gt.69.0 .and. E_emerge.lt.91.0 )then
				// std::cout << i,Erange_det[i],(range_det-DetThick); 
				// endif 
			} 

			// calculate the energy of the emergant ejectile
			det_E = linear(Erange_det[i-1],Erange_det[i],E_det[i-1],E_det[i],(range_det-DetThick)); 
			det_dE = E_emerge - det_E; 

			// **************************************************
			// Angular Straggling

			if(HitGeo==1){
				//strag_dE_cyl(A,Z,DetThick,DetPhi_cyl(detHitCyl),DetRad_cyl[DetHitCyl],theta_Emerge,phi_Emerge,E_emerge,theta_E,phi_E,DetRadL,conv_det); 
			}
			else if(HitGeo==2){ 
			 	strag_dE_plan(A,Z,DetThick,DetZ_plan[DetHitPlan],theta_Emerge,phi_Emerge,E_emerge,theta_E,phi_E,DetRadL,conv_det); 
			}
			else if(HitGeo==3){ 
				//theta_E = theta_Emerge; 
				//phi_E = phi_Emerge; 
		
				// Put in Angular Straggling for the annular case. Only needed to shift strips in dE and E case
				// Always have dE strip number anyway, so not so important.
			} 
		} // Particle stops in dE check

		//------------------------------------------------------------------; 
		// E DETECTOR
		// Calculate the ejectile's range in the E detector material
		if(det_E>0){// SDP 1.5 MeV threshold
			i = 0; 
			while(E_det[i] < det_E){ i = i + 1; } 

			range_det = linear(E_det[i-1],E_det[i],Erange_det[i-1],Erange_det[i],det_E); 

			// Calculate the thickness the ejectile sees in the E layer
			//  Note that this uses the effective straggling angle rather than the real angle!!!
			//  Correct this by getting AngStrag to return the real angle as well!!

			if(HitGeo==1){ 
				//det_thick_cyl(theta_E,phi_E,DetPhi_cyl(detHitCyl),Ethick_cyl(detHitCyl),DetThick); 
			}
			else if(HitGeo==2){ 
				det_thick_plan(theta_E,phi_E,Ethick_plan[DetHitPlan],DetThick); 
			} 
			else if(HitGeo==3){ 
				//det_thick_plan(theta_E,phi_E,Ethick_ann(detHitAnn),DetThick); 
			}

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

		if(HitGeo == 1){
			// Cylindrical Geometry
			/*resolution_cyl(det_dE,theta_Emerge,phi_Emerge,beamspot,targ_angle,Eres_dE_cyl(detHitCyl),Pres_dE_cyl(detHitCyl),
				       Pres_dE_en_cyl(detHitCyl),detLength_cyl[DetHitCyl],DetRad_cyl[DetHitCyl],det_dE_final,theta_dE_final); 

			/*if(det_E==0.0){
				det_E_final = 0.0; 
				theta_E_final = 0.0; 
			} 
			else{ 
				resolution_cyl(det_E,theta_E,phi_E,beamspot,targ_angle,Eres_E_cyl(detHitCyl),Pres_E_cyl(detHitCyl),
					       Pres_E_en_cyl[detHitCyl],detLength_cyl[DetHitCyl],DetRad_cyl[DetHitCyl],det_E_final,theta_E_final); 
			}*/
		}
		else if(HitGeo==2){ 
			// Planar Geometry
			det_dE_final = det_dE; 
			theta_dE_final = theta_Emerge; 
			det_E_final = det_E; 
			theta_E_final = theta_Emerge; 

			resolution_plan(det_dE,theta_Emerge,phi_Emerge,Eres_dE_plan[DetHitPlan],Pres_dE_plan[DetHitPlan],Pres_dE_en_plan[DetHitPlan],
					DetZ_plan[DetHitPlan],xPlane[DetHitPlan],yPlane[DetHitPlan],DetXMin[DetHitPlan][StripHitPlan],DetYMin[DetHitPlan][StripHitPlan],
					detWidth_plan[DetHitPlan],detLength_plan[DetHitPlan],Nstrips_plan[DetHitPlan],StripHitPlan,det_dE_final,theta_dE_final,phi_dE_final); 
				
			if(det_E==0.0){
				det_E_final = 0.0; 
				theta_E_final = 0.0; 
			} 
			else{ 
				resolution_plan(det_E,theta_E,phi_E,Eres_E_plan[DetHitPlan],Pres_E_plan[DetHitPlan],Pres_E_en_plan[DetHitPlan],
						DetZ_plan[DetHitPlan],xPlane[DetHitPlan],yPlane[DetHitPlan],DetXMin[DetHitPlan][StripHitPlan],DetYMin[DetHitPlan][StripHitPlan],
						detWidth_plan[DetHitPlan],detLength_plan[DetHitPlan],Nstrips_plan[DetHitPlan],StripHitPlan,det_E_final,theta_E_final,phi_E_final); 
			} 
		}
		else if(HitGeo==3){ 
			// Annular Geometry
			/*det_dE_final = det_dE; 
			theta_dE_final = theta_Emerge; 
			det_E_final = det_E; 
			theta_E_final = theta_Emerge; 

			resolution_ann(det_dE,theta_Emerge,phi_Emerge,Eres_dE_ann(detHitAnn),DetZ_Ann(DetHitAnn),DetInner_ann(detHitAnn),DetOuter_ann(detHitAnn),
				       DetPhi_ann(detHitAnn),Nstrips_Ann(detHitAnn),StripHitAnn,det_dE_final,theta_dE_final,phi_dE_final); 

			if(det_E==0.0){
				det_E_final = 0.0; 
				theta_E_final = 0.0; 
			} 
			else{ 
				resolution_ann(det_E,theta_E,phi_E,Eres_E_ann(detHitAnn),DetZ_Ann(DetHitAnn),DetInner_ann(detHitAnn),DetOuter_ann(detHitAnn),
					       DetPhi_ann(detHitAnn),Nstrips_Ann(detHitAnn),StripHitAnn,det_E_final,theta_E_final,phi_E_final); 
			} */
		} // Over HitGeometry
	} // HitCounter = 1 
} 

#endif
