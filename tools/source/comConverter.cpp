
#include <fstream>

#include "vandmc_core.hpp"

#include "comConverter.hpp"

comConverter::comConverter() : length(0) { }

comConverter::comConverter(const char *fname) : length(0) { 
	load(fname);
}

bool comConverter::load(const char *fname){
	com.clear();
	ejectLab.clear();
	recoilLab.clear();
	length = 0;

	std::ifstream file(fname);
	if(!file.good()) return false;

	double comVal, ejectLabVal, recoilLabVal;
	while(true){
		file >> comVal >> ejectLabVal >> recoilLabVal;
		if(file.eof()) break;
		com.push_back(comVal);
		ejectLab.push_back(ejectLabVal);
		recoilLab.push_back(recoilLabVal);
	}

	length = com.size();
	return (length != 0);
}

double comConverter::convertEject2lab(const double &com_){
	double retval = -9999;
	Interpolate(com_, retval, com.data(), ejectLab.data(), length);
	return retval;
}

double comConverter::convertRecoil2lab(const double &com_){
	double retval = -9999;
	Interpolate(com_, retval, com.data(), recoilLab.data(), length);
	return retval;
}
