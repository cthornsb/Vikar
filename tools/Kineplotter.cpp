// Cory Thornsberry
// 05/11/2014
// Kineplotter

#define PLOTTER_VERSION "05.11.2014A"

#include "/home/cory/Research/ROOTlib/FSExecutable.hpp"
#include "/home/cory/Research/ROOTlib/FSPlotter.hpp"

using namespace FightingSquirrel;

int main(int argc, char* argv[]){	
	std::cout << "------------------------------------------------------------\n";
	std::cout << "-                      Kineplotter                         -\n";
	std::cout << "------------------------------------------------------------\n";

	if(!(argc >= 3)){
		std::cout << " Error: Insufficienct number of parameters given, aborting\n";
		return 1;
	}
	
	if(strcmp(argv[1],"version") == 0){ // Display version information
		std::cout << " KinePlotter v" << PLOTTER_VERSION << " by Cory Thornsberry\n";
		std::cout << "  FightingSquirrel Library Dependencies\n";
		std::cout << "   Framework\tv" << FSVersions::Framework() << std::endl;
		std::cout << "   Filler\tv" << FSVersions::Filler() << std::endl;
		return 0;
	}

	char* chicken[0]; TApplication* rootapp = new TApplication("rootapp",0,chicken);

	// Main
	Filler filler;
	unsigned short num_states = short(atol(argv[2]));
	unsigned int count = filler.Initialize(argv[1], num_states);
	if(count > 0){
		filler.GetGraph()->GetXaxis()->SetTitle("Lab Theta (deg)");
		filler.GetGraph()->GetYaxis()->SetTitle("Ejectile Energy (MeV)");
		filler.GetGraph()->SetTitle("Ejectile Energy vs. Lab Angle");
		filler.Draw("AP");
		filler.Pause();
	}
	else{
		std::cout << " Note: Found no entries\n";
	}

	return 0;
}
