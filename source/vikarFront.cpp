// vikar_front.cpp
// Front-End for Vikar 3.12
// Cory Thornsberry
// Feb. 15th, 2014

#include <string>
#include <fstream>
#include <iostream>

int main(){
	unsigned short version = 312;
	unsigned short idummy, Nstates, Nelements;
	double ddummy; 
	std::string input;

	std::cout << "\n####       #### ######## ####     ###         ###       #########\n"; 
	std::cout << " ##         ##     ##     ##     ##          ## ##       ##     ##\n";
	std::cout << "  ##       ##      ##     ##    ##          ##   ##      ##     ##\n";
	std::cout << "  ##       ##      ##     ##   ##           ##   ##      ##   ##\n";
	std::cout << "   ##     ##       ##     #####            #########     #####\n";
	std::cout << "   ##     ##       ##     ##   ##          ##     ##     ##   ##\n";
	std::cout << "    ##   ##        ##     ##    ##        ##       ##    ##    ##\n";
	std::cout << "    ##   ##        ##     ##     ##       ##       ##    ##     ##\n";
	std::cout << "     ## ##         ##     ##      ##     ##         ##   ##      ##\n";
	std::cout << "      ###       ######## ####      ###  ####       #### ####      ###\n";

	std::cout << "\n VIKAR 3.12\n"; 
	std::cout << " ==  ==  ==  ==  == \n\n"; 
	
	std::cout << " Welcome to VIKAR << the Virtual Instrumentation\n"; 
	std::cout << " for Kinematics And Reactions program front-end program\n\n"; 

	std::cout << " How about a nice cup of tea?\n"; 
	
	// Read in the reaction type to be simulated
	std::cout << " No? Well then, let's get started shall we?\n";

	std::cout << "\n Enter the filename for your input file: "; 
	std::cin >> input;
	std::ofstream file91(input.c_str());
	
	if(!file91.good()){
		std::cout << " Error: Problem opening the output file\n";
		return 1;
	}

	file91 << "1\n"; // Silly but kept in for backwards compatibility
	file91 << version/100.0 << "\n"; 

	// Read in the reaction type to be simulated
	std::cout << "\n Is your reaction to be 2 or 3 body?\n"; 
	std::cout << " (NB - only 2 body reactions supported presently): "; 
	std::cin >> idummy;

	while(idummy != 2){ // .and.body.ne.3)
		std::cout << "\n Please select a 2 or 3 body reaction: "; 
		std::cin >> idummy;
	}  
	
	file91 << idummy << std::endl;

	//------------------------------------------------------------------------
	// Read in the reaction data
	std::cout << "\n Reaction Data\n"; 
	std::cout << " -------------\n"; 
	
	std::cout << "\n Enter your beam species' Z: "; 
	std::cin >> input; file91 << input << std::endl;
	std::cout << " Enter your beam species' A: ";
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n Enter your target nuclide's Z and A.\n"; 
	std::cout << " Note that these are the details for\n"; 
	std::cout << " the reaction and kinematics only. You\n"; 
	std::cout << " will subsequently be required to enter\n"; 
	std::cout << " details of the complete target\n"; 
	std::cout << " composition << for energy loss and\n"; 
	std::cout << " straggling calcualtions.\n"; 
	
	std::cout << "\n Enter your target nuclide's Z: "; 
	std::cin >> input; file91 << input << std::endl;
	std::cout << " Enter your target nuclide's A: ";
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n Enter your ejectile's Z: "; 
	std::cin >> input; file91 << input << std::endl;
	std::cout << " Enter your ejectile's A: ";
	std::cin >> input; file91 << input << std::endl;

	// Read in the beam energy
	std::cout << "\n Enter the beam energy (MeV): "; 
	std::cin >> input; file91 << input << std::endl;

	// Read in the beam spot size
	std::cout << "\n Enter the beam spot size FWHM (mm): "; 
	std::cin >> input; file91 << input << std::endl;

	// Read in the beam energy spread
	std::cout << "\n Enter the beam enery spread FWHM (MeV): "; 
	std::cin >> input; file91 << input << std::endl;

	// Read in data about the states populated in the final nucleus
	std::cout << "\n Enter the Q-value for the reaction where all\n"; 
	std::cout << " particles are in their ground state (in MeV): "; 
	std::cin >> input; file91 << input << std::endl;

	// Read in data about the states populated in the final nucleus
	std::cout << "\n Enter the number of states populated\n"; 
	std::cout << " in the residual: "; 
	std::cin >> Nstates; 
	file91 << Nstates << std::endl;

	std::cout << "\n Enter the excitation energy of the states\n"; 
	std::cout << " populated in the residual << in order of\n"; 
	std::cout << " increasing excitation energy (in MeV).\n"; 

	for(unsigned short i = 0; i < Nstates; i++){
		std::cout << "  State " << i+1 << ": "; 
		std::cin >> input; file91 << input << std::endl;
	} 

	std::cout << "\n Do you wish to supply individual angular\n"; 
	std::cout << " distributions for ejectiles for the states populated\n"; 
	std::cout << " in the residual << or do you wish to assume isotropy\n"; 
	std::cout << " in the CoM? (0 = isotropic reactions << 1 = supply distributions): "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = isotropic reactions << 1 = supply distributions): "; 
		std::cin >> idummy;
	} 

	file91 << idummy << std::endl;

	if(idummy == 1){
		for(unsigned short i = 0; i < Nstates; i++){
			std::cout << "  State " << i+1 << ": "; 
			std::cin >> input; file91 << input << std::endl;
		} // over Nstates
	} 

	// End of Reaction Data section

	//------------------------------------------------------------------------
	//
	// Start of Energy Loss section
	//
	//------------------------------------------------------------------------
 
	std::cout << "\n --------------------------------------------------\n"; 
	std::cout << " Range Tables\n"; 
	std::cout << " --------------------------------------------------\n"; 
	std::cout << "\n The energy loss of the beam particles << ejectiles\n"; 
	std::cout << " and recoils can be calculated (independently)\n"; 
	std::cout << " from external range tables or using internal routines.\n"; 
	std::cout << " Note that the use of the internal routines is recommended.\n"; 
	std::cout << " Note also that currently only SRIM output is\n"; 
	std::cout << " supported for external range tables.\n"; 

	// Read in the target details
	std::cout << "\n Target Properties\n"; 
	std::cout << " -----------------\n"; 
	std::cout << "\n Enter details of youe target composition << for\n"; 
	std::cout << " the purposes of energy loss and straggling\n"; 
	std::cout << " calcualtions. Currently << only single layer targets\n"; 
	std::cout << " are supported. Enter the number of elements in your target molecule: "; 
	std::cin >> Nelements;
	file91 << Nelements << std::endl;

	for(unsigned short i = 0; i < Nelements; i++){
		std::cout << "  Z of Element " << i+1 << ": "; 
		std::cin >> input; file91 << input << std::endl;
		std::cout << "  A of Element " << i+1 << ": "; 
		std::cin >> input; file91 << input << std::endl;
		std::cout << "  No. of Element " << i+1 << " per molecule: "; 
		std::cin >> input; file91 << input << std::endl;
	} 

	// Read in target thickness
	std::cout << "\n Enter the target thickness (mg/cm^2): "; 
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n The target density is required for internal range\n"; 
	std::cout << " tables (this variable will not be used if you later\n";
	std::cout << " select external tables).\n"; 
	std::cout << " Some common densities are given below.\n"; 
	std::cout << " CH2 = 0.93 g/cm^3\n"; 
	std::cout << " CD2 = 1.063 g/cm^3\n"; 
	std::cout << " Mylar = 1.39 g/cm^3\n"; 
	std::cout << " LiO = 2.01 g/cm^3\n"; 
	std::cout << " C = 2.265 g/cm^3\n"; 
	std::cout << " Al = 2.70 g/cm^3\n"; 
	std::cout << " Mg = 1.70 g/cm^3\n"; 
	std::cout << " Si = 2.33 g/cm^3\n"; 
	std::cout << " Pb = 11.35 g/cm^3\n"; 
	std::cout << " Au = 19.32 g/cm^3\n"; 
	std::cout << " Enter the target density (g/cm^3): "; 
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n Enter the target angle (degrees).\n"; 
	std::cout << " A target angle of 0 corresponds to\n"; 
	std::cout << " a target plane perpendicular to the\n"; 
	std::cout << " beam axis. Rotation is clockwise\n"; 
	std::cout << " about the y axis << as viewed from above: "; 
	std::cin >> input; file91 << input << std::endl;

	//------------------------------------------------------------------------
	// Code for processing range tables for BEAM PARTICLES in TARGET
	// Currently only SRIM output is supported for external tables

	// Determine the energy loss method for beam particles
	std::cout << "\n Beam particle Tables\n"; 
	std::cout << " --------------------\n"; 
	std::cout << "\n Enter (1) to read in your own range table or\n"; 
	std::cout << " (0) to rely on the internal calculations: "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = internal range tables << 1 = external range tables): "; 
		std::cin >> idummy;
	} 

	file91 << idummy << std::endl;

	// If external range tables are to be used for beam particles
	if(idummy == 1){
		std::cout << "\n Enter the file name (and path if necessary)\n";
		std::cout << " of the SRIM file to be read: "; 
		std::cin >> input; file91 << input << std::endl;
	}

	std::cout << "\n An excellent choice, if I may be so bold. Now,\n";
	std::cout << "\n Ejectiles in Target Tables\n"; 
	std::cout << " --------------------------\n"; 
	std::cout << "\n Enter (1) to read in your own range table or\n"; 
	std::cout << " (0) to rely on the internal calculations: "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = internal range tables << 1 = external range tables): "; 
		std::cin >> idummy;
	} 
	
	file91 << idummy << std::endl;

	// If external range tables are to be used for ejectiles
	if(idummy == 1){
		std::cout << "\n Enter the file name (and path if necessary)";
		std::cout << "\n of the SRIM file for ejectiles in TARGET: "; 
		std::cin >> input; file91 << input << std::endl;
	} // Ejectile tables.eq.1

	//------------------------------------------------------------------------
	// Code for processing range tables for RECOILS in TARGET
	// Currently only SRIM output is supported for external tables

	std::cout << "\n Recoils in Target Tables\n"; 
	std::cout << " ------------------------\n"; 
	std::cout << "\n Enter (1) to read in your own range table or\n"; 
	std::cout << " (0) to rely on the internal calculations: "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = internal range tables << 1 = external range tables): "; 
		std::cin >> idummy;
	} 
	
	file91 << idummy << std::endl;

	// If external range tables are to be used for recoils
	if(idummy == 1){
		// Read in for all particles in the target material
		std::cout << "\n Enter the file name (and path if necessary)";
		std::cout << " of the SRIM file for recoils in TARGET: "; 
		std::cin >> input; file91 << input << std::endl;
	} // Recoil tables.eq.1

	//------------------------------------------------------------------------
	// Code for processing range tables for EJECTILES in DETECTOR
	// Currently only SRIM output is supported for external tables

	std::cout << "\n Ejectiles in Detector Tables\n"; 
	std::cout << " ----------------------------\n"; 
	std::cout << "\n Currently all detectors are made of Si\n"; 
	std::cout << "\n Enter (1) to read in your own range table or\n"; 
	std::cout << " (0) to rely on the internal calculations: "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = internal range tables << 1 = external range tables): "; 
		std::cin >> idummy;
	}  
	
	file91 << idummy << std::endl;
	
	// If external range tables are to be used for ejectiles
	if(idummy == 1){
		std::cout << "\n Enter the file name (and path if necessary)\n";
		std::cout << " of the SRIM file for ejectiles in DETECTOR: "; 
		std::cin >> input; file91 << input << std::endl;
	} // Ejectile tables.eq.1

	//------------------------------------------------------------------------
	// Code for processing range tables for RECOILS in DETECTOR
	// Currently only SRIM output is supported for external tables

	std::cout << "\n Recoils in Detector Tables\n"; 
	std::cout << " --------------------------\n"; 
	std::cout << "\n Currently all detectors are made of Si\n"; 
	std::cout << "\n Enter (1) to read in your own range table or\n"; 
	std::cout << " (0) to rely on the internal calculations: "; 
	std::cin >> idummy;

	while(idummy != 0 && idummy != 1){
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = internal range tables << 1 = external range tables): "; 
		std::cin >> idummy;
	} 
	
	file91 << idummy << std::endl;

	// If external range tables are to be used for recoils
	if(idummy == 1){
		std::cout << "\n Enter the file name (and path if necessary)\n";
		std::cout << " of the SRIM file for recoils in DETECTOR: "; 
		std::cin >> input; file91 << input << std::endl;
	} // recoil tables.eq.1

	//------------------------------------------------------------------------
	//
	// End of Energy Loss section
	//
	//------------------------------------------------------------------------

	// Read in detector locations, distances, resolutions, and particle type sensitivities.
	std::cout << "\n --------------------------------------------------\n"; 
	std::cout << " Detection System\n"; 
	std::cout << " --------------------------------------------------\n"; 
	std::cout << "\n The simulation can be run to simulate a float\n"; 
	std::cout << " detection system << with detectors with finite\n"; 
	std::cout << " coverage and resolution << or a perfect detection\n"; 
	std::cout << " system can be used << in which all particles are\n"; 
	std::cout << " detected << with perfect resolution.\n"; 
	std::cout << "\n Would you like to include the detection setup in the simulation?\n"; 
	std::cout << " Enter (1) for float detectors << or (0) for a perfect detector: "; 
	std::cin >> idummy; 

	while(idummy != 0 && idummy != 1){ 
		std::cout << "\nUser input incorrect. Please select\n"; 
		std::cout << "(0 = perfect detector << 1 = float detector): "; 
		std::cin >> idummy;
	} 

	file91 << idummy << std::endl;

	// If 'real' detectors are being used
	if(idummy == 1){
		unsigned short DetSetup = 0;
		std::cout << "\n A floatist << are you?\n"; 
		std::cout << "\n If you have a simple detection system << you can\n"; 
		std::cout << " enter it here. Alternatively << you can specify the\n"; 
		std::cout << " path and filename of a detector setup file.\n"; 
		std::cout << " (Note that currently only file based setup is accepted)."; 
		std::cout << "\n Enter (0) to manually specify << or (1) to use a setup file: "; 
		std::cin >> DetSetup;

		while(DetSetup != 0 && DetSetup != 1){
			std::cout << "\nUser input incorrect. Please select\n"; 
			std::cout << "(0 = manual setup << 1 = setup file): "; 
			std::cin >> DetSetup;
		} 
		
		file91 << DetSetup << std::endl;

		// If the detector setup is being read in from an external file
		if(DetSetup == 1){
			std::cout << "\n Enter the file name (and path if necessary)\n";
			std::cout << " of the detector setup file for cylindrical\n"; 
			std::cout << " geometry: "; 
			std::cin >> input; file91 << input << std::endl;

			std::cout << "\n Enter the file name (and path if necessary)\n";
			std::cout << " of the detector setup file for planar geometry: "; 
			std::cin >> input; file91 << input << std::endl;

			std::cout << "\n\n Enter the file name (and path if necessary)\n"; 
			std::cout << " of the detector setup file for annular geometry: "; 
			std::cin >> input; file91 << input << std::endl;
		} // DetSetup
	} // PerfectDet

	// End of detector set up section

	//------------------------------------------------------------------------
	//
	// Beam Quality Determination
	//
	//------------------------------------------------------------------------

	// Read in here the beam spot size, the beam divergence, and the energy spread

	// Put in something here that allows you to choose whether the source of events
	// is restricted to the origin, or whether to move the detectors around in space
	// for each event

	//------------------------------------------------------------------------
	// Read in the number of detected events wanted
	std::cout << "\n Statistics\n"; 
	std::cout << " ----------\n"; 
	std::cout << "\n Enter the desired number of detected events: "; 
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n Experiment Sensitivity\n"; 
	std::cout << " ----------------------\n"; 
	std::cout << "\n The particles to be detected are specified here\n"; 
	std::cout << " in order to determine what constitutes a hit << \n"; 
	std::cout << " and to avoid the full simulation of particles not\n"; 
	std::cout << " being detected. Are ejectiles detected? (Enter 1 for yes << 0 for no): "; 
	std::cin >> idummy; file91 << idummy << std::endl;

	if(idummy == 1){ // Det_eject = true
		std::cout << "\n Are recoils detected? (Enter 1 for yes << 0 for no): "; 
		std::cin >> idummy; file91 << idummy << std::endl;
		if(idummy == 1){ // Det_eject = Det_recoil = true
			std::cout << "\n Is a coincidence between ejectiles and recoils\n"; 
			std::cout << " required << or are only singles required?\n"; 
			std::cout << " (Enter 1 for yes << 0 for no): "; 
			std::cin >> input; file91 << input << std::endl;
		}
	}
	else{ // Det_eject = false
		std::cout << "\n Are recoils detected? (Enter 1 for yes << 0 for no): "; 
		std::cin >> idummy; file91 << idummy << std::endl;
	}

	std::cout << "\n Write out undetected ejectiles?\n"; 
	std::cout << " (For low detection efficiencies << answering\n"; 
	std::cout << " yes will slow the simulation appreciably)\n"; 
	std::cout << " Enter 1 for yes << 0 for no: "; 
	std::cin >> input; file91 << input << std::endl;

	std::cout << "\n Write out undetected recoils?\n"; 
	std::cout << " (For low detection efficiencies << answering\n"; 
	std::cout << " yes will slow the simulation appreciably)\n"; 
	std::cout << " Enter 1 for yes << 0 for no: "; 
	std::cin >> input; file91 << input << std::endl;

	file91.close();
	return 0;
}
