#include "structures.h"

///////////////////////////////////////////////////////////
// RecoilObject
///////////////////////////////////////////////////////////
void RecoilObject::Append(const double &hitX_, const double &hitY_, const double &hitZ_, const double &hitTheta_, const double &hitPhi_,
						  const double &qdc_, const double &tof_, const double &faceX_, const double &faceY_, const double &faceZ_, const int &loc_){
	recoil_hitX.push_back(hitX_);
	recoil_hitY.push_back(hitY_);
	recoil_hitZ.push_back(hitZ_);
	recoil_hitTheta.push_back(hitTheta_);
	recoil_hitPhi.push_back(hitPhi_);
    recoil_qdc.push_back(qdc_);
    recoil_tof.push_back(tof_);
    recoil_faceX.push_back(faceX_);
    recoil_faceY.push_back(faceY_);
    recoil_faceZ.push_back(faceZ_);
    recoil_loc.push_back(loc_);
    recoil_mult++;
}

void RecoilObject::Zero(){
	if(recoil_mult == 0){ return ; } // Structure is already empty
	recoil_hitX.clear();
	recoil_hitY.clear();
	recoil_hitZ.clear();
	recoil_hitTheta.clear();
	recoil_hitPhi.clear();
    recoil_qdc.clear();
    recoil_tof.clear();
    recoil_faceX.clear();
    recoil_faceY.clear();
    recoil_faceZ.clear();
    recoil_loc.clear();
	recoil_mult = 0;
}

///////////////////////////////////////////////////////////
// EjectObject
///////////////////////////////////////////////////////////
void EjectObject::Append(const double &hitX_, const double &hitY_, const double &hitZ_, const double &hitTheta_, const double &hitPhi_,
						 const double &qdc_, const double &tof_, const double &faceX_, const double &faceY_, const double &faceZ_, const int &loc_){
	eject_hitX.push_back(hitX_);
	eject_hitY.push_back(hitY_);
	eject_hitZ.push_back(hitZ_);
	eject_hitTheta.push_back(hitTheta_);
	eject_hitPhi.push_back(hitPhi_);
    eject_qdc.push_back(qdc_);
    eject_tof.push_back(tof_);
    eject_faceX.push_back(faceX_);
    eject_faceY.push_back(faceY_);
    eject_faceZ.push_back(faceZ_);
    eject_loc.push_back(loc_);
    eject_mult++;
}

void EjectObject::Zero(){
	if(eject_mult == 0){ return ; } // Structure is already empty
	eject_hitX.clear();
	eject_hitY.clear();
	eject_hitZ.clear();
	eject_hitTheta.clear();
	eject_hitPhi.clear();
    eject_qdc.clear();
    eject_tof.clear();
    eject_faceX.clear();
    eject_faceY.clear();
    eject_faceZ.clear();
    eject_loc.clear();
	eject_mult = 0;
}
