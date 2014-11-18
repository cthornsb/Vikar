#include <iostream>
#include <fstream>

#include "vikar_core.h"
#include "materials.h"
#include "detectors.h"
#include "structures.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

double momentum(double energy_, double mass_){
	return (1.0/(3E8))*std::sqrt(energy_*energy_ + 2.0*mass_*energy_);
}

int main(int argc, char *argv[]){
	gSystem->Load("libTree");
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

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
	
	/*TFile *file = new TFile("VIKAR.root","READ");
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

	rootapp->Delete();*/

	/*TFile *file = new TFile("VIKAR.root","READ");
	TTree *tree = (TTree*)file->Get("VIKAR");	

	TCanvas *can = new TCanvas("can");
	can->cd();

	std::vector<unsigned int> loc;
	std::vector<double> qdc;
	TBranch *qdc_b, *loc_b;
	
	tree->SetMakeClass(1);
	tree->SetBranchAddress("recoil_qdc", &qdc, &qdc_b);
	tree->SetBranchAddress("recoil_loc", &loc, &loc_b);
	
	TH2D *hist = new TH2D("hist", "Phoswich Fast QDC vs. Slow", 100, 0, 28, 100, 0, 28);
	hist->SetStats(false);
	hist->GetXaxis()->SetTitle("E (MeV)");
	hist->GetYaxis()->SetTitle("dE (MeV)");
	
	std::vector<unsigned int>::iterator iterl;
	std::vector<double>::iterator iterq;
	
	double qdc1, qdc2;
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		qdc1 = 0.0; qdc2 = 0.0;
		for(iterl = loc.begin(), iterq = qdc.begin(); iterl != loc.end() && iterq != qdc.end(); iterl++, iterq++){
			if(*iterl == 42){ // Fast
				qdc1 = *iterq;
			}
			else if(*iterl == 43){ // Slow
				qdc2 = *iterq;
			}
		}
		
		if(qdc1 != 0.0 && qdc2 != 0.0){
			hist->Fill(qdc2, qdc1);
		}
	}

	hist->Draw("COLZ");
	can->WaitPrimitive();

	file->Close();
	can->Close();*/

	TFile *file = new TFile("Elastic.root","READ");
	TTree *tree = (TTree*)file->Get("VIKAR");	

	TCanvas *can = new TCanvas("can");
	can->cd();

	double mass7Be = 4*(proton_RME) + 3*(neutron_RME); // MeV
	double massD = proton_RME + neutron_RME; // MeV

	std::vector<double> recoilE, recoilX, recoilY, recoilZ;
	std::vector<double> ejectE, ejectX, ejectY, ejectZ;
	std::vector<double> reactE, reactX, reactY, reactZ;
	TBranch *recoilE_b, *recoilX_b, *recoilY_b, *recoilZ_b;
	TBranch *ejectE_b, *ejectX_b, *ejectY_b, *ejectZ_b;
	TBranch *reactE_b, *reactX_b, *reactY_b, *reactZ_b;
	
	tree->SetMakeClass(1);
	tree->SetBranchAddress("eject_qdc", &ejectE, &ejectE_b);
	tree->SetBranchAddress("eject_hitX", &ejectX, &ejectX_b);
	tree->SetBranchAddress("eject_hitY", &ejectY, &ejectY_b);
	tree->SetBranchAddress("eject_hitZ", &ejectZ, &ejectZ_b);
		
	tree->SetBranchAddress("recoil_qdc", &recoilE, &recoilE_b);
	tree->SetBranchAddress("recoil_hitX", &recoilX, &recoilX_b);
	tree->SetBranchAddress("recoil_hitY", &recoilY, &recoilY_b);
	tree->SetBranchAddress("recoil_hitZ", &recoilZ, &recoilZ_b);
	
	tree->SetBranchAddress("reactE", &reactE, &reactE_b);
	tree->SetBranchAddress("trajectoryX", &reactX, &reactX_b);
	tree->SetBranchAddress("trajectoryY", &reactY, &reactY_b);
	tree->SetBranchAddress("trajectoryZ", &reactZ, &reactZ_b);
	
	TH2D *hist = new TH2D("hist", "Beam Lab Angle Reconstruction", 100, 0.0, 5.0, 100, 0.0, 1.0);
	hist->SetStats(false);
	hist->GetYaxis()->SetTitle("|Pbeam - (Peject+Precoil)|/|Pbeam|");
	hist->GetXaxis()->SetTitle("θabs (deg)");	

	std::vector<double>::iterator iterEE, iterEX, iterEY, iterEZ;
	std::vector<double>::iterator iterRE, iterRX, iterRY, iterRZ;
	std::vector<double>::iterator iterReactE, iterReactX, iterReactY, iterReactZ;
	double maximum_theta = -9999;
	double beam_theta;
	Vector3 eject_momentum;
	Vector3 recoil_momentum;
	Vector3 beam_momentum;
	Vector3 sum_vector;
	
	std::cout << " Processing " << tree->GetEntries() << " events\n";
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		iterEE = ejectE.begin(); iterEX = ejectX.begin(); iterEY = ejectY.begin(); iterEZ = ejectZ.begin();
		iterRE = recoilE.begin(); iterRX = recoilX.begin(); iterRY = recoilY.begin(); iterRZ = recoilZ.begin();
		iterReactE = reactE.begin(); iterReactX = reactX.begin(); iterReactY = reactY.begin(); iterReactZ = reactZ.begin();
		
		// Construct the momentum vectors
		// For the ejectile...
		Vector3 ejectile(*iterEX, *iterEY, *iterEZ);
		Cart2Sphere(ejectile, eject_momentum);
		eject_momentum.axis[0] = momentum(*iterEE, massD);
		Sphere2Cart(eject_momentum);
		
		// For the recoil...
		Vector3 recoil(*iterRX, *iterRY, *iterRZ);
		Cart2Sphere(recoil, recoil_momentum);
		recoil_momentum.axis[0] = momentum(*iterRE, mass7Be);
		Sphere2Cart(recoil_momentum);
		
		// And for the beam particle...
		Vector3 beam(*iterReactX, *iterReactY, *iterReactZ);
		Cart2Sphere(beam, beam_momentum);
		beam_momentum.axis[0] = momentum(*iterReactE, mass7Be);
		beam_theta = beam_momentum.axis[1];
		Sphere2Cart(beam_momentum);
		
		sum_vector = eject_momentum + recoil_momentum;
		if(beam_theta > maximum_theta){ maximum_theta = beam_theta; }		
		
		hist->Fill(beam_theta*rad2deg, (beam_momentum-sum_vector).Length()/beam_momentum.Length());
	}
	
	std::cout << " Done! Maximum theta = " << maximum_theta*rad2deg << std::endl;

	hist->Draw("COLZ");
	can->WaitPrimitive();

	file->Close();
	can->Close();

	rootapp->Delete();	

	return 0;
}
