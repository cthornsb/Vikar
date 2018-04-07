#ifndef DATA_PACK_HPP
#define DATA_PACK_HPP

#include "Structures.h"

class TFile;
class TTree;

class dataPack{
  public:
	TFile *file;
	TTree *tree;
	bool init;

	double offsetX, offsetY, offsetZ;
	double trajX, trajY, trajZ;
	double Ereact, Eeject, Erecoil;
	double labTheta, labPhi;
	double comAngle;
	int labBin;

	MonteCarloStructure MCARLOdata;

	dataPack();
	
	dataPack(std::string fname_, unsigned int Nbins_, bool write_rxn_=false);

	~dataPack();

	bool IsInit();

	bool Open(std::string fname_, bool write_rxn_);

	bool Close();
	
	void Zero();
};

#endif
