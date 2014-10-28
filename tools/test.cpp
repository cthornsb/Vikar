#include <iostream>
#include <fstream>

#include "vikar_core.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

int main(){
	gSystem->Load("libTree");
	//char* dummy[0]; 
	//TApplication* rootapp = new TApplication("rootapp",0,dummy);
	
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

	TFile *file = new TFile("VIKAR.root","READ");
	TTree *tree = (TTree*)file->Get("VIKAR");

	TCanvas *can = new TCanvas("can");
	can->cd();
		
	for(unsigned int i = 0; i < 42; i++){
		std::cout << " Processing eject_loc = " << i << std::endl;
		std::stringstream stream; stream << i;
		tree->Draw("recoil_faceY:recoil_faceX>>(100,-0.025,0.015,100,-0.02,0.02)",("eject_loc=="+stream.str()).c_str(),"COLZ");
		can->Print("out.gif+10");
	}

	can->Close();

	//rootapp->Delete();

	return 0;
}
