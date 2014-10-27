// structures.h
// Cory Thornsberry

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "TObject.h"

#include <vector>

///////////////////////////////////////////////////////////
// RecoilObject
///////////////////////////////////////////////////////////
class RecoilObject : public TObject {
  public:
    std::vector<double> recoil_hitX, recoil_hitY, recoil_hitZ, recoil_hitTheta, recoil_hitPhi;
    std::vector<double> recoil_qdc, recoil_tof, recoil_faceX, recoil_faceY, recoil_faceZ;
    std::vector<int> recoil_loc;
    unsigned int recoil_mult;

    RecoilObject(){ recoil_mult = 0; }
    
    // Add an entry to the data vector
    // Calling this method will mark the event as valid
    void Append(const double &hitX_, const double &hitY_, const double &hitZ_, const double &hitTheta_, const double &hitPhi_,
    			const double &qdc_, const double &tof_, const double &faceX_, const double &faceY_, const double &faceZ_, const int &loc_);
    
    // Zero the data structure
    void Zero();
    
    ClassDefNV(RecoilObject, 1); // Recoil
};

///////////////////////////////////////////////////////////
// EjectObject
///////////////////////////////////////////////////////////
class EjectObject : public TObject {
  public:
    std::vector<double> eject_hitX, eject_hitY, eject_hitZ, eject_hitTheta, eject_hitPhi;
    std::vector<double> eject_qdc, eject_tof, eject_faceX, eject_faceY, eject_faceZ;
    std::vector<int> eject_loc;
    unsigned int eject_mult;

	EjectObject(){ eject_mult = 0; }
    
    // Add an entry to the data vector
    // Calling this method will mark the event as valid
    void Append(const double &hitX_, const double &hitY_, const double &hitZ_, const double &hitTheta_, const double &hitPhi_, 
    			const double &qdc_, const double &tof_, const double &faceX_, const double &faceY_, const double &faceZ_, const int &loc_);
    
    // Zero the data structure
    void Zero();
    
    ClassDefNV(EjectObject, 1); // Eject
};

///////////////////////////////////////////////////////////
// ReactionObject
///////////////////////////////////////////////////////////
class ReactionObject : public TObject {
  public:
  	std::vector<double> reactE;
    std::vector<double> reactX, reactY, reactZ;
    std::vector<double> trajectoryX, trajectoryY, trajectoryZ;
    unsigned int reaction_mult;

	ReactionObject(){ reaction_mult = 0; }
    
    // Add an entry to the data vector
    // Calling this method will mark the event as valid
    void Append(const double &reactE_, const double &reactX_, const double &reactY_, const double &reactZ_,
    			const double &trajectoryX_, const double &trajectoryY_, const double &trajectoryZ_);
    
    // Zero the data structure
    void Zero();
    
    ClassDefNV(ReactionObject, 1); // Reaction
};

#endif
