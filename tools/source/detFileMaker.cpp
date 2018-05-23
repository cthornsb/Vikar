


#include <iostream>
#include <fstream>
#include <string>

#include "vandmc_core.hpp"

int main(){
	double radius;
	double spacing;
	double startAngle;
	int Ndet;

	std::cout << " Enter detector radius (m): "; std::cin >> radius;
	std::cout << " Enter detector spacing (deg): "; std::cin >> spacing;
	std::cout << " Enter number of detectors: "; std::cin >> Ndet;
	std::cout << " Enter angle of 1st det (deg): "; std::cin >> startAngle;

	std::string filename;
	std::cout << " Enter detector filename: "; std::cin >> filename;

	std::ofstream file(filename.c_str());
	if(!file.good()){
		std::cout << "  ERROR! Failed to open output file \"" << filename << "\"!\n";
		return -1;
	}

	file << "#x\ty\tz\ttheta\tphi\tpsi\ttype\tsybtype\tlength\twidth\tdepth\tmaterial\n";

	double x, y, z, angle;
	for(int i = 0; i < Ndet; i++){
		angle = startAngle + spacing*i;
		Sphere2Cart(radius, angle*deg2rad, 0, x, y, z);
		file << x << "\t" << y << "\t" << z << "\t" << angle*deg2rad << "\t0\t0\tvandle\tsmall\n";
	}

	std::cout << "  Done! Wrote file \"" << filename << "\".\n";

	file.close();

	return 0;
}
