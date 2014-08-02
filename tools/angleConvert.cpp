// angleConvert.cpp
// Cory Thornsberry
// April 25, 2014
// Convert an xmgrace x-section file from CoM angle to lab angle

#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>

#include "vikar_core.h"

// Return true if substr is found in the character array
bool FindSubstring(char *input, char* substr, unsigned short len){
	bool found = false;
	unsigned short i = 0, j;
	while(input[i] != '\0'){
		if(input[i] == substr[0]){ 
			found = true; j = 0;
			while(j < len){
				if(input[i+j] == '\0'){ return false; }
				if(input[i+j] != substr[j]){ 
					found = false;
					break; 
				}
				j++;
			}
			
			if(found){ return true; }
		}
		i++;
	}
	return false;
}

// Find a single character in the character array
bool FindChar(char *input, char substr){
	for(unsigned short i = 0; i < 256; i++){
		if(input[i] == substr){ return true; }
	}
	return false;
}

// Parse a string from an xmgrace file
// return true if values are found and false otherwise
bool ParseString(char *input, float &theta, float &sigma){
	if(FindChar(input,'@') || FindChar(input,'#')){ return false; }
	else{
		// No xmgrace junk found, this must be data
		bool value = false;
		bool xvalue = true;
		std::string value_str = "";
		for(unsigned short i = 0; i < 256; i++){
			if(input[i] != ' '){ value = true; }
			if(value){
				if(input[i] != ' '){ value_str += input[i]; }
				else{
					value = false;
					if(xvalue){ 
						theta = atof(value_str.c_str()); 
						value_str = "";
						xvalue = false;
					}
					else{ sigma = atof(value_str.c_str()); }
				}
			}
		}
	}
}

void help(){
	std::cout << " Help-------\n";
	std::cout << "  Syntax: angleConvert [filename]\n";
	std::cout << "   filename must be the path to the xmgrace (fort) file\n\n";
	
	std::cout << " Program will ask the user a series of questions including\n";
	std::cout << "  information about the kinematics reaction which is taking\n";
	std::cout << "  place in the system.\n\n";
	
	std::cout << " The number of states corresponds to the number of allowed\n";
	std::cout << "  excitation states of the recoil. If the program encounters\n";
	std::cout << "  more or less data in the fort file than the number of states\n";
	std::cout << "  provided, it will print a warning message, but will not crash.\n\n";
	
	std::cout << " Providing fewer states than in the input file will cause the\n";
	std::cout << "  extra data in the file to be ignored since the recoil excitation\n";
	std::cout << "  is not provided. Providing more states will effectively do nothing.\n\n";
}

int main(int argc, char* argv[]){
	std::ifstream inFile;
	if(argc >= 2){
		// Load input file
		if(strcmp(argv[1],"help") == 0){
			help();
			return 0;
		}
		
		inFile.open(argv[1]);
		if(!inFile.good()){
			std::cout << " Error: Failed to load input file " << argv[1] << std::endl;
			return 1;
		}
	}
	else{
		std::cout << " Error: No input filename given\n";
		std::cout << "  Try 'angleConvert help' for correct syntax\n";
		return 1;
	}
	
	Kindeux kind;	
	double Ebeam, Mbeam, Mtarg, Meject, Q;
	unsigned short num_states;
	double *RecoilEx;
	std::cout << "  Enter beam Energy (MeV): "; std::cin >> Ebeam;
	std::cout << "  Enter beam Mass (A): "; std::cin >> Mbeam;
	std::cout << "  Enter target Mass (A): "; std::cin >> Mtarg;
	std::cout << "  Enter g.s. Q-value (MeV): "; std::cin >> Q;
	std::cout << "  Enter ejectile Mass (A): "; std::cin >> Meject;
	std::cout << "  Enter number of states: "; std::cin >> num_states;
	
	if(num_states > 0){
		RecoilEx = new double[num_states];
		for(unsigned short i = 0; i < num_states; i++){
			std::cout << "   Enter recoil excitation for state " << i+1 << " (MeV): "; 
			std::cin >> RecoilEx[i];
		}
	}
	else{
		std::cout << " Error: Expected non-zero number of excited states\n";
		return 1;
	}
	
	std::cout << std::endl;
	kind.Initialize(Mbeam, Mtarg, ((Mbeam+Mtarg)-Meject), Meject, Q);
	
	std::ofstream outFile("angular.out");
	double temp;
	float theta, sigma;
	char line[256];
	char end[] = "END";
	int count = 0, badcount = 0;
	unsigned short state = 0;
	while(true){
		// Load values from input file
		inFile.getline(line,256);
		if(inFile.eof()){ break; }
		
		// Encountered more data than expected
		if(state+1 > num_states){
			std::cout << " Warning! Found more than " << num_states << " state(s) in file\n\n";
			break;
		}
			
		// Parse string and look for data
		if(!FindSubstring(line,end,3)){
			if(ParseString(line, theta, sigma)){
				// String contains data and not xmgrace junk
				temp = WrapValue(kind.ConvertAngle2CoM(Ebeam, 0.0, RecoilEx[state], theta*deg2rad)*rad2deg, 0.0, 360.0);				
				if(temp >= 0.0){ outFile << temp << "\t" << sigma << "\n"; }
				else{ 
					std::cout << " Cannot convert theta = " << theta << " sigma = " << sigma << std::endl;
					badcount++; 
				}
				count++;
			}
			else{
				// String contains xmgrace info
				outFile << line << "\n";
			}
		}
		else{ 
			state++;
			outFile << line << "\n"; 
		}
	}
	inFile.close();
	outFile.close();
	
	// Encountered less data than expected
	if(state < num_states){
		std::cout << " Warning! Received " << num_states-state << " extraneous state(s)\n";
		std::cout << "  Found only " << state << " states in file\n\n";
	}
	
	std::cout << " Converted " << count << " data points\n";
	std::cout << " Lost " << badcount << " bad data points\n";
	std::cout << " Wrote output file angular.out\n\n";
	
	delete[] RecoilEx;
	
	return 0;
}
