#include <iostream>
#include <fstream>

#include "vikar_core.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

int main(){
	gSystem->Load("libTree");
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	
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

	/*TFile *file = new TFile("VIKAR.root","READ");
	TTree *tree = (TTree*)file->Get("VIKAR");

	TCanvas *can = new TCanvas("can");
	can->cd();
		
	for(unsigned int i = 0; i < 42; i++){
		std::cout << " Processing eject_loc = " << i << std::endl;
		std::stringstream stream; stream << i;
		tree->Draw("recoil_faceY:recoil_faceX>>(100,-0.025,0.015,100,-0.02,0.02)",("eject_loc=="+stream.str()).c_str(),"COLZ");
		can->Print("out.gif+10");
	}

	can->Close();*/
	
	TFile *file = new TFile("VIKAR.root","READ");
	TTree *tree = (TTree*)file->Get("VIKAR");	

	TCanvas *can = new TCanvas("can");
	can->cd();

	std::vector<unsigned int> loc;
	std::vector<double> x, y, z;
	TBranch *x_b, *y_b, *z_b, *loc_b;
	
	tree->SetMakeClass(1);
	tree->SetBranchAddress("recoil_hitX", &x, &x_b);
	tree->SetBranchAddress("recoil_hitY", &y, &y_b);
	tree->SetBranchAddress("recoil_hitZ", &z, &z_b);
	tree->SetBranchAddress("recoil_loc", &loc, &loc_b);
	
	TH1D *hist1 = new TH1D("hist1", "Theta (deg)", 180, 0.0, 20.0);
	TH1D *hist2 = new TH1D("hist2", "Phi (deg)", 360, 0.0, 360.0);
	
	std::vector<unsigned int>::iterator iterl;
	std::vector<double>::iterator iterx, itery, iterz;
	double r, theta, phi;
	double max_theta = -999;
	double max_phi = -999;
	
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		for(iterx = x.begin(), itery = y.begin(), iterz = z.begin(), iterl = loc.begin();
		  iterx != x.end() && itery != y.end() && iterz != z.end() && iterl != loc.end(); iterx++, itery++, iterz++, iterl++){
			if(*iterl == 42){
				Cart2Sphere(*iterx, *itery, *iterz, r, theta, phi);
				theta *= rad2deg; phi *= rad2deg;
				
				hist1->Fill(theta);
				hist2->Fill(phi);
				
				if(theta > max_theta){ max_theta = theta; }
				if(phi > max_phi){ max_phi = phi; }
			}
		}
	}
	
	std::cout << " Maximum theta: " << max_theta << " deg\n";
	std::cout << " Maximum phi: " << max_phi << " deg\n";

	can->Divide(2);
	can->cd(1);
	hist1->Draw();
	can->cd(2);
	hist2->Draw();
	can->WaitPrimitive();

	rootapp->Delete();

	return 0;
}
