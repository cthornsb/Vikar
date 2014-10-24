#include <iostream>
#include <fstream>

#include "vikar_core.h"

int main(){
	Vector3 vector(0.0, 0.0, 0.5); // Start along the z-axis
	Matrix3 matrix;
	
	std::ofstream output("test.dat");
	
	double theta_step = pi/10.0;
	double phi_step = 2*pi/10.0;
	double current_theta;
	double current_phi;
	for(unsigned int i = 0; i <= 10; i++){
		current_theta = i*theta_step;
		for(unsigned int j = 0; j <= 10; j++){
			vector = Vector3(0.0, 0.0, 0.5);
			current_phi = j*phi_step;
			
			matrix.SetRotationMatrixSphere(current_theta, current_phi);	
			matrix.Transform(vector);
			output << vector.axis[0] << "\t" << vector.axis[1] << "\t" << vector.axis[2] << "\n";
		}
	}
	
	std::cout << " Done!\n";
	output.close();

	return 0;
}
