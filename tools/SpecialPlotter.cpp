// Cory Thornsberry
// 08/19/2013
// Plotter

#define PLOTTER_VERSION "04.09.2014A"

#include <cmath>
#include "/home/cory/Research/ROOTlib/FSExecutable.hpp"
#include "/home/cory/Research/ROOTlib/FSExtractor.hpp"
#include "/home/cory/Research/ROOTlib/FSPlotter.hpp"

using namespace FightingSquirrel;

// Return a random number between low and high
double frand(double low, double high){
	return low+(double(rand())/RAND_MAX)*(high-low);
}

// Mimic the fortran rand() function
double frand(){
	return double(rand())/RAND_MAX;
}

double rndgauss0(double w){
	// rndgauss0, a cut down version of rndgauss1; 
	// returns a random number with FWHM w centred at 0;
	
	double t, tsq;  
	const double c0=2.515517, c1=0.802853, c2=0.010328; 
	const double d1=1.432788, d2=0.189269, d3=0.001308; 
	const double widthfact=0.424628450; 

	if(w == 0.0){ return 0.0; } 

	t = frand(); 

	if(t > 0.5){ t = t-0.5; }
	if(t < 1e-30){ t = 11.46380587; }
	else{ 
		tsq = -log(t*t); 
		t = std::sqrt(tsq); 
	} 
	
	//     compute inverse by equn 26.2.23 
	t=t-(c0+c1*t+c2*tsq)/(1.00+d1*t+(d2+d3*t)*tsq);
	
	//     now randomize x in positive and negative direction
	if(frand() > 0.5){ t = -t; }
	
	//     now correct for standard deviation
	return widthfact*w*t; 
} 

void help(const char* option=nullptr){
	std::cout << "HELP:\n";	
	std::cout << " SpecialPlotter {filename.root} {0-7}\n";
	std::cout << "  Histogram Types-----------------------------\n";
	std::cout << "  0 : (VIKARoutput) Plot Lab Theta vs. Energy for the Ejectile\n";
	std::cout << "  1 : (VIKARoutput) Plot Lab Theta vs. ToF for the Ejectile\n";
	std::cout << "  2 : (VIKARoutput) Plot Lab Theta vs. Lab Phi for the Ejectile\n";
	std::cout << "  3 : (VIKARoutput) Plot Bar# vs. Hit_Y for the Ejectile\n\n";
	std::cout << "  4 : (RECOILoutput) Plot Lab Theta vs. Energy for the Recoil\n";
	std::cout << "  5 : (RECOILoutput) Plot Lab Theta vs. Lab Phi for the Recoil\n\n";
	std::cout << "  6 : (DEBUGoutput) Plot Lab Theta vs. Energy for debug data\n";
	std::cout << "  7 : (DEBUGoutput) Plot Lab Theta vs. Lab Phi for debug data\n\n";
}

int main(int argc, char* argv[]){	
	std::cout << "------------------------------------------------------------\n";
	std::cout << "-                     SpecialPlotter                       -\n";
	std::cout << "------------------------------------------------------------\n";

	if(!(argc >= 2)){
		std::cout << " Error: No parameters given, aborting\n";
		std::cout << "  Try: SpecialPlotter [help], for more information\n";
		return 1;
	}
	
	if(strcmp(argv[1],"help") == 0){
		if(argc >= 3){ help(argv[2]); }
		else{ help(); }
		return 0;
	}
	else if(strcmp(argv[1],"version") == 0){ // Display version information
		std::cout << " SpecialPlotter v" << PLOTTER_VERSION << " by Cory Thornsberry\n";
		std::cout << "  FightingSquirrel Library Dependencies\n";
		std::cout << "   Framework\tv" << FSVersions::Framework() << std::endl;
		std::cout << "   Extractor\tv" << FSVersions::Extractor() << std::endl;
		std::cout << "   Plotter2D\tv" << FSVersions::Plotter2D() << std::endl;
		return 0;
	}
	else if(!(argc >= 3)){
		std::cout << " Error: Insufficient number of parameters provided, aborting\n";
		std::cout << "  Try: SpecialPlotter [help], for more information\n";
		return 1;
	}

	char* chicken[0]; TApplication* rootapp = new TApplication("rootapp",0,chicken);
	bool theta_energy = false; 
	bool theta_tof = false;
	bool energy_bar = false;
	bool tof_bar = false;
	bool theta_phi = false;
	bool bar_y = false;
	bool recoil_theta_energy = false;
	bool recoil_theta_phi = false;
	bool debug_theta_energy = false;
	bool debug_theta_phi = false;

	Extractor data;
	if(!data.LoadFile(argv[1])){
		std::cout << " Error: The input file failed to open, aborting\n";
		return 45;
	}
	TTree *tree_ = data.Extract<TTree>("Tree");
	if(tree_ == nullptr){
		std::cout << " Error: Failed to load the input tree, aborting\n";
		return 46;
	}

	if(strcmp(argv[2],"0") == 0){ theta_energy = true; } // Plot Lab Theta vs. Energy
	else if(strcmp(argv[2],"1") == 0){ theta_tof = true; } // Plot Lab Theta vs. ToF
	else if(strcmp(argv[2],"2") == 0){ energy_bar = true; } // Plot Energy vs. Bar#
	else if(strcmp(argv[2],"3") == 0){ tof_bar = true; } // Plot ToF vs. Bar#	
	else if(strcmp(argv[2],"4") == 0){ theta_phi = true; } // Plot Lab Theta vs. Lab Phi
	else if(strcmp(argv[2],"5") == 0){ bar_y = true; } // Plot Bar# vs. Hit_Y
	else if(strcmp(argv[2],"6") == 0){ recoil_theta_energy = true; } // Plot Lab Theta vs. Energy for the Recoil
	else if(strcmp(argv[2],"7") == 0){ recoil_theta_phi = true; } // Plot Lab Theta vs. Lab Phi for the Recoil
	else if(strcmp(argv[2],"8") == 0){ debug_theta_energy = true; } // Plot Lab Theta vs. Energy for debug data
	else if(strcmp(argv[2],"9") == 0){ debug_theta_phi = true; } // Plot Lab Theta vs. Lab Phi for debug data
	else{
		std::cout << " Error: Unrecognized option '" << argv[2] << "', aborting\n";
		std::cout << "  Try: SpecialPlotter [help], for more information\n";
		return 3;
	}

	// Main
	Plotter2D Plot;
	unsigned int count;
	if(theta_energy){ 
		Plot.Initialize(tree_->GetName(),tree_,180,twoTuple<double>(0.0,180.0),400,twoTuple<double>(0.0,4.0)); 
		std::cout << " Found " << Plot.FillHist("Col03","Col05") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (deg)");	
		Plot.GetHist()->GetYaxis()->SetTitle("Neutron Energy (MeV)");	
	}
	else if(theta_tof){ 
		Plot.Initialize(tree_->GetName(),tree_,180,twoTuple<double>(0.0,180.0),400,twoTuple<double>(0.0,150));
		std::cout << " Found " << Plot.FillHist("Col03","Col06") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Neutron ToF (ns)");
		Plot.GetHist()->GetYaxis()->SetTitle("Lab Theta (deg)");
	}
	else if(energy_bar){ 
		Plot.Initialize(tree_->GetName(),tree_,400,twoTuple<double>(0.0,4.0),42,twoTuple<double>(0.0,42.0));
		std::cout << " Found " << Plot.FillHist("Col05","Col07") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Neutron Energy (MeV)");
		Plot.GetHist()->GetYaxis()->SetTitle("Bar#");
	}
	else if(tof_bar){ 
		Plot.Initialize(tree_->GetName(),tree_,333,twoTuple<double>(0.0,150.0),42,twoTuple<double>(0.0,42.0));
		std::cout << " Found " << Plot.FillHist("Col06","Col07") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Neutron ToF (ns)");
		Plot.GetHist()->GetYaxis()->SetTitle("Bar#");
		
		/*// Temporary gamma flash
		double time, bar;
		for(unsigned short i = 0; i < 10000; i++){
			time = frand(1.6667, 1.9433); // Random gamma
			bar = frand(0, 42); // Random bar
			time += rndgauss0(3.0); // Smear by 3 ns
			Plot.GetHist()->Fill(time,bar);
		}*/
	}
	else if(theta_phi){ 
		Plot.Initialize(tree_->GetName(),tree_,180,twoTuple<double>(0.0,180.0),360,twoTuple<double>(0.0,360.0)); 
		std::cout << " Found " << Plot.FillHist("Col03","Col04") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (deg)");
		Plot.GetHist()->GetYaxis()->SetTitle("Lab Phi (deg)");
	}
	else if(bar_y){ 
		Plot.Initialize(tree_->GetName(),tree_,42,twoTuple<double>(0.0,42.0),21,twoTuple<double>(-1.0,1.0)); 
		std::cout << " Found " << Plot.FillHist("Col07","Col10") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Bar#");
		Plot.GetHist()->GetYaxis()->SetTitle("Y-Hit (m)");
	}
	else if(recoil_theta_energy){ 
		Plot.Initialize(tree_->GetName(),tree_,200,twoTuple<double>(0.0,0.2),400,twoTuple<double>(10.0,25.0)); 
		std::cout << " Found " << Plot.FillHist("Col00","Col02") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (rad)");
		Plot.GetHist()->GetYaxis()->SetTitle("Recoil Energy (MeV)");
	}
	else if(recoil_theta_phi){ 
		Plot.Initialize(tree_->GetName(),tree_,200,twoTuple<double>(0.0,0.2),400,twoTuple<double>(0.0,6.28318)); 
		std::cout << " Found " << Plot.FillHist("Col00","Col01") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (rad)");
		Plot.GetHist()->GetYaxis()->SetTitle("Lab Phi (rad)");
	}
	else if(debug_theta_energy){ 
		Plot.Initialize(tree_->GetName(),tree_,400,twoTuple<double>(0.0,3.14159),400,twoTuple<double>(0.0,12.0));
		std::cout << " Found " << Plot.FillHist("Col00","Col02") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (rad)");
		Plot.GetHist()->GetYaxis()->SetTitle("Neutron Energy (MeV)");
	}
	else if(debug_theta_phi){ 
		Plot.Initialize(tree_->GetName(),tree_,200,twoTuple<double>(0.0,3.14159),400,twoTuple<double>(0.0,6.28318));
		std::cout << " Found " << Plot.FillHist("Col00","Col01") << " entries\n";
		Plot.GetHist()->GetXaxis()->SetTitle("Lab Theta (rad)");
		Plot.GetHist()->GetYaxis()->SetTitle("Lab Phi (rad)");
	}
	
	// Fill a 2-d histogram with values
	Plot.GetHist()->GetYaxis()->SetTitleOffset(1.2);
	Plot.MakeDrawable();
	Plot.Draw("COLZ",false);
	Plot.Pause();

	return 0;
}
