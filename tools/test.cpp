#include <iostream>
#include <fstream>

#include "vikar_core.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1D.h"

int main(){
	gSystem->Load("libTree");
	/*Vector3 vector(0.0, 0.0, 0.5); // Start along the z-axis
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
	output.close();*/

	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

	AngularDist dist;
	dist.Initialize("/home/cory/Research/VANDLE/Vikar/angdist.dat", 0.0, 0.0, 0.0);
	std::cout << " Num Points: " << dist.GetNumPoints() << std::endl;
	std::cout << " X-section: " << dist.GetReactionXsection() << " mb\n";
	
	TH1D *hist = new TH1D("hist", "hist", 180, 0, 180.0);
	for(unsigned int i = 0; i < 100000; i++){
		hist->Fill(dist.Sample());
	}

	TCanvas *can = new TCanvas("can");
	can->cd();
	hist->Draw();
	can->WaitPrimitive();
	can->Close();

	rootapp->Delete();

	return 0;
}
