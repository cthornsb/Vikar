// PulseAnalyzer.cpp
// C. Thornsberry
// Aug. 25th, 2014
// Load detector pulses and analyze them using various integration methods
// SYNTAX: Viewer {root_filename} {root_treename} {root_branch} {method#} [debug]

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TApplication.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <cmath>

const double PIXIE_TIME_RES = 4.0; // In ns
const double ROOT2PI = std::sqrt(2*3.14159);

struct peak{
	unsigned int left, max, min, cfd;
	double max_value, min_value;
	
	peak(){
		left = 0; max = 0; min = 0; cfd = 0; 
		max_value = -9999.0; min_value = 9999.0;
	}
	
	peak(unsigned int left_, unsigned int max_, unsigned int min_, unsigned int cfd_, double max_value_, double min_value_){
		left = left_; max = max_; min = min_; cfd = cfd_; 
		max_value = max_value_; min_value = min_value_;
	}
	
	std::string print(){
		std::stringstream output;
		output << "left = " << left << ", cfd = " << cfd;
		output << ", max_bin = " << max << ", max_value = " << max_value;
		output << ", min_bin = " << min << ", min_value = " << min_value; 
		return output.str();
	}
};

double interpolate(double x, double x1, double y1, double x2, double y2){
	if(x == x1){ return y1; }
	if(x == x2){ return y2; }
	return (y1 + ((y2-y1)/(x2-x1))*(x-x1)); 
}

unsigned int preprocess(double *input, unsigned int start, unsigned int stop, unsigned int size, peak *first, peak *second, bool &valley, bool debug=false){
	// Do baseline correction. Use the range to the left of the pulse.
	// The right side of the pulse will have a slow component
	double base_sum = 0.0;
	unsigned int base_max = start + (unsigned int)(0.1*size); // First 10% of range
	for(unsigned int i = start; i < base_max; i++){
		base_sum += input[i];
	}

	// Correct the baseline
	double base_line = base_sum/(base_max-start);
	for(unsigned int i = start; i < stop; i++){
		input[i] = input[i] - base_line;
	}

	// Find the global maximum
	double maximum = -9999.0;
	unsigned int maximum_bin = 0;
	for(unsigned int i = start; i < stop; i++){
		if(input[i] > maximum){ 
			maximum = input[i];
			maximum_bin = i;
		}
	}

	// Peak-finder (there should be more than 1)
	unsigned int leading_edge, peak_bin;
	double back_sub_left = 0.0;
	double back_sub_right = 0.0;
	double slope;
	bool up_slope = false;
	bool down_slope = false;
	bool found_peak = false;
	unsigned short num_peaks = 0;
	for(unsigned int i = start; i < stop-1; i++){
		// Subtract background (10% of maximum)
		back_sub_left = input[i] - 0.1*maximum;
		back_sub_right = input[i+1] - 0.1*maximum;
		if(back_sub_left < 0.0){ back_sub_left = 0.0; }
		if(back_sub_right < 0.0){ back_sub_right = 0.0; }
		
		// Calculate the slope
		slope = back_sub_right - back_sub_left;
		if(!up_slope && slope > 0.0){ // Moving to up_slope marks leading edge of pulse
			up_slope = true; 
			down_slope = false;
			leading_edge = i;
		}
		else if(up_slope && slope < 0.0){ // Moving from up_slope to down_slope marks a pulse peak
			down_slope = true; 
			up_slope = false;
			found_peak = true;
			num_peaks++;
			peak_bin = i;
		}
		else if(down_slope && slope > 0.0){ // Moving from down_slope to non-negative slope marks trailing edge (reset peak-finder)
			up_slope = false;
			down_slope = false;
		}
		
		if(found_peak){
			found_peak = false;
			if(num_peaks == 1){
				first->left = leading_edge;
				first->max = peak_bin;
				first->max_value = input[peak_bin];
			}
			else if(num_peaks == 2){
				second->left = leading_edge;
				second->max = peak_bin;
				second->max_value = input[peak_bin];
				break; // That's enough!
			}
		}
	}

	// Find the valley if it exists
	if(num_peaks > 1){ // Fast and slow component, valley will lie between the first and second maxima
		for(unsigned int i = first->max; i < second->max; i++){
			if(input[i] < first->min_value){
				first->min_value = input[i];
				first->min = i;
			}
		}
		valley = true;
	}
	else{ // Only the fast component. The valley does not exist
		valley = false;
	}
	
	// Print debug information
	if(debug){ 
		std::cout << "Global: baseline = " << base_line << ", maximum_bin = " << maximum_bin;
		std::cout << ", maximum = " << maximum << ", num_peaks = " << num_peaks << std::endl;
		std::cout << "Peak: " << first->print() << std::endl;
	}
	
	return num_peaks;
}

// Integrate a pulse and return the short and long integrals
unsigned int integrate(std::vector<int> &arr, unsigned int pulse_size, std::vector<double> &s_int, std::vector<double> &l_int, unsigned short method, bool debug, TCanvas *canvas1){
	if(debug && (!canvas1)){ debug = false; }
	unsigned int size = arr.size();
	s_int.clear(); l_int.clear();
	if(size == 0){ return 0; }
	if(size % pulse_size != 0){ return 0; }
	unsigned int num_pulses = size/pulse_size; // Expects size to be divisible by pulse_size	
	
	// Copy the pulse into a new array so the original pulse is unmodified
	double *darr = new double[size];
	unsigned int count = 0;
	for(std::vector<int>::iterator iter = arr.begin(); iter != arr.end(); iter++){
		if(count >= size){ break; }
		darr[count] = (double)(*iter);
		count++;
	}
	
	// Variables for TGraph
	double *x1 = new double[pulse_size];
	double *x2 = new double[pulse_size];
	const double *x_val1 = new double[pulse_size];
	const double *y_val1 = new double[pulse_size];
	const double *x_val2 = new double[pulse_size];
	const double *y_val2 = new double[pulse_size];
	
	TGraph *graph1 = new TGraph(pulse_size, x_val1, y_val1);
	TGraph *graph2 = new TGraph(pulse_size, x_val2, y_val2);
	for(unsigned int i = 0; i < pulse_size; i++){ x1[i] = i; x2[i] = i; }
		
	if(debug){
		graph1->SetLineColor(2); graph1->SetLineWidth(2);
		graph2->SetLineColor(4); graph2->SetLineWidth(2);
		canvas1->cd();
	}
	
	bool has_valley;
	unsigned int peak_count;
	unsigned int start, stop;
	peak first_peak;
	peak second_peak;
	double s, l;
	
	//if(num_pulses != 2){ return 0; }
	//for(unsigned int pulse = 1; pulse < num_pulses; pulse++){
	for(unsigned int pulse = 0; pulse < num_pulses; pulse++){
		start = pulse*pulse_size;
		stop = (pulse+1)*pulse_size;

		peak_count = preprocess(darr, start, stop, pulse_size, &first_peak, &second_peak, has_valley, debug);

		if(debug || method == 5 || method == 6){
			for(unsigned int i = start; i < stop; i++){ 
				graph1->SetPoint(i-start, x1[i-start], darr[i]); 
			}
		}

		if(peak_count == 0){ // Failed to find a peak
			if(debug){ std::cout << "No peaks found, skipping entry\n"; }
			continue;
		}
		if((method >= 4) && peak_count > 1){ // Multiple peaks will throw off the fitting functions
			if(debug){ std::cout << "Multiple peaks, skipping entry\n"; }
			continue;
		}
	
		// Do short and long integrals
		s = 0.0; l = 0.0;			
		if(method == 0){ // Slope inversion
			// To the right of the maximum, we expect a negative slope
			// Stop short integral when the slope goes non-negative
			bool calc_short = true;
			for(unsigned int i = first_peak.cfd; i < first_peak.min; i++){
				// Integrate using trapezoid rule
				if(calc_short){ // Waiting on positive slope to stop calculating short integral
					if(i < first_peak.max || darr[i+1]-darr[i] <= 0.0){ s += (darr[i] + darr[i+1])/2.0; } // Negative slope, still on side of pulse
					else{ calc_short = false; }
				}
				l += (darr[i] + darr[i+1])/2.0;
			}
		} // method == 0
		else if(method == 1){ // Fixed pulse height
			// Stop short integral when we reach 10% of the pulse height
			// Find the 10% of maximum point
			unsigned int ten_percent_bin = 0;
			for(unsigned int i = first_peak.max; i < stop; i++){
				if(darr[i] <= 0.1*first_peak.max_value){ 
					ten_percent_bin = i;
					break;
				}
			}
			
			if(debug){ std::cout << " method = 1, 10p_bin = " << ten_percent_bin << std::endl; }	
			
			for(unsigned int i = first_peak.cfd; i < first_peak.min; i++){
				// Integrate using trapezoid rule with arbitrary bin width of 1
				if(i < ten_percent_bin){ s += (darr[i] + darr[i+1])/2.0; } // Stop calculating short below 10% maximum
				l += (darr[i] + darr[i+1])/2.0;
			}
		} // method == 1
		else if(method == 2){ // Search for decay-side slope variation
			unsigned int sstop = 0;
			unsigned int minimum_slope_bin;
			double minimum_slope;
			if(has_valley){ // A valley was found. This is a clear separation between fast and slow
				sstop = first_peak.min;
				if(debug){
					std::cout << " method = 2, sstop = " << sstop << std::endl;			
				}
			}
			else{ // No valley was found. Search for slope variation to the right of the fast pulse		
				double darr_prime[pulse_size];
				for(unsigned int i = start; i < stop-1; i++){
					darr_prime[i-start] = darr[i+1] - darr[i];
				}
				darr_prime[pulse_size-1] = 0.0;

				// Find the most negative slope, this will be the characteristic slope of the pulse
				double minimum_slope = 9999.0;
				unsigned int minimum_slope_bin = 0;
				for(unsigned int i = start; i < stop; i++){
					if(darr_prime[i] < minimum_slope){
						minimum_slope = darr_prime[i-start];
						minimum_slope_bin = i-start;
					}
				}

				// Find the bin in which to stop short integration
				for(unsigned int i = minimum_slope_bin; i < stop; i++){
					if(darr_prime[i] >= 0.5*minimum_slope){
						sstop = i;
						break;
					}
				}
				
				if(debug){
					std::cout << " method = 2, slope_bin = " << minimum_slope_bin << ", min_slope = " << minimum_slope;
					std::cout << ", sstop = " << sstop << std::endl;			
					for(unsigned int i = start; i < stop; i++){ // 1st derivative
						graph2->SetPoint(i-start, x2[i-start], darr_prime[i-start]);
					}
					canvas1->Clear();
					graph2->Draw();
					graph1->Draw("SAME");
					canvas1->Update();
					canvas1->WaitPrimitive();
				}	
			}
			
			// Integrate using trapezoid rule
			for(unsigned int i = start; i < stop; i++){
				if(i < sstop){ s += (darr[i] + darr[i+1])/2.0; }
				l += (darr[i] + darr[i+1])/2.0;
			}
			
			if(debug && s == 0){ 
				std::cout << "peak_count = " << peak_count << std::endl;
				std::cout << "slope_bin = " << minimum_slope_bin << ", min_slope = " << minimum_slope;
				std::cout << ", sstop = " << sstop << ", lstop = " << first_peak.min << std::endl;
				std::cout << first_peak.print() << std::endl << std::endl;
			}
		} // method == 2

		if(debug){ std::cout << " s = " << s << ", l = " << l << std::endl; }
		s_int.push_back(s);
		l_int.push_back(l);
	} // Over num_pulses
	
	graph1->Delete(); 
	graph2->Delete();
	delete[] darr; delete[] x1; delete[] x2; delete[] x_val1; 
	delete[] y_val1; delete[] x_val2; delete[] y_val2;
	return num_pulses;
}

// For compilation
int main(int argc, char* argv[]){
	std::string method_names[3] = {"slope inversion (S vs. L)", "fixed height (S vs. L)", "decay-side slope variation (S vs. L)"};

	if(argc < 4){
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./analyzer {filename} {treename} {method#} [debug]\n";
		for(unsigned short i = 0; i < 3; i++){
			std::cout << "  Method " << i << ": " << method_names[i] << std::endl;
		}
		return 1;
	}
		
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	unsigned short method = atol(argv[3]);
	std::cout << " Using analysis method " << method << ": " << method_names[method] << std::endl;
	if(method > 2){
		std::cout << " Encountered undefined method (" << method << "), aborting\n";
		return 1;
	}
	
	bool debug = false;
	if(argc > 4){
		if(strcmp(argv[4], "debug") == 0){ 
			debug = true; 
			std::cout << " DEBUGGING...\n";
		}
	}

	// Branch variables
	std::vector<int> wave;
	int type;
	unsigned int wave_size = 0;

	TFile *file = new TFile(argv[1], "READ");
	if(file->IsZombie()){
		std::cout << " Failed to load the input file '" << argv[1] << "'\n";
		return 1;
	}
	TTree *tree = (TTree*)file->Get(argv[2]);
	if(!tree){
		std::cout << " Failed to load the input tree '" << argv[2] << "'\n";
		file->Close();
		return 1;
	}
	tree->SetMakeClass(1);

	TBranch *b_wave, *b_type;
	tree->SetBranchAddress("pulse", &wave, &b_wave);
	tree->SetBranchAddress("type", &type, &b_type);
	
	if(!b_wave){
		std::cout << " Failed to load the input branch '" << argv[3] << "'\n";
		file->Close();
		return 1;
	}

	double short_value, long_value;
	TFile *out_file = NULL;
	TTree *out_tree = NULL;
	if(!debug){
		std::cout << " Opening output file 'analyzer.root'\n";
		out_file = new TFile("analyzer.root", "RECREATE");
		
		// Standard tree
		out_tree = new TTree(argv[2], "Pulse analysis tree");
		out_tree->Branch("short", &short_value);
		out_tree->Branch("long", &long_value);
		out_tree->Branch("type", &type);
	}

	// Get the pulse size
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if(wave.size() == 0){ continue; }
		else{ 
			wave_size = wave.size();
			break; 
		}
	}
	std::cout << " Using wave size " << wave_size << std::endl;

	// Canvas
	TCanvas *can1 = new TCanvas("can1", "canvas");
	can1->cd();

	// Histogram
	TH2D *hist = NULL;
	unsigned int status_update = 100000;
	if(method == 0){ // slope inversion
		hist = new TH2D("hist", method_names[0].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 1){ // constant percentage
		hist = new TH2D("hist", method_names[1].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 2){ // slope variation
		hist = new TH2D("hist", method_names[2].c_str(), 50, 55, 85, 50, 0, 1);
		hist->GetXaxis()->SetTitle("L/1000 (arb. units)");
		hist->GetYaxis()->SetTitle("(L-S)/L (arb. units)");
	}
	hist->SetStats(false);

	std::vector<double> short_integral, long_integral;
	std::vector<double>::iterator iter1, iter2;
	unsigned int count = 0;
	
	long time_holder1 = 0;
	double time_holder2 = 0.0;
	clock_t cpu_time = clock();
	time_t real_time;
	time(&real_time);
	
	unsigned int num_entries = tree->GetEntries();
	if(debug){ num_entries = 100; }
	
	std::cout << " Processing " << num_entries << " entries\n";
	for(unsigned int i = 0; i < num_entries; i++){
		tree->GetEntry(i);
		if(i % status_update == 0){ 
			time_holder1 += (long)(clock()-cpu_time);
			time_holder2 += difftime(time(NULL), real_time);
			cpu_time = clock();
			time(&real_time);
			if(i != 0){ 
				if(time_holder2 < 1){ time_holder2 = 1; } // Prevent zero time remaining
				std::cout << " Entry no. " << i << ", CPU = " << ((float)time_holder1)/CLOCKS_PER_SEC << " s, REAL = " << time_holder2;
				std::cout << " s, REMAIN = " << ((float)(num_entries-i))*(time_holder2/i) << " s\n";
			}
		}
		count += integrate(wave, wave_size, short_integral, long_integral, method, debug, can1);
		for(iter1 = short_integral.begin(), iter2 = long_integral.begin(); iter1 != short_integral.end() && iter2 != long_integral.end(); iter1++, iter2++){
			short_value = (*iter1);
			long_value = (*iter2);
			hist->Fill(long_value/1000.0, 1.0-(short_value/long_value));
			if(!debug){ out_tree->Fill(); }
		}
	}
	
	std::cout << " Found " << count << " pulses in " << num_entries << " tree entries\n";
	if(!debug){ 
		std::cout << " Wrote " << out_tree->GetEntries() << " entries to the output tree\n";
	}

	can1->cd();
	hist->Draw("COLZ");
	can1->Update();
	can1->WaitPrimitive();

	if(!debug){
		out_file->cd();
		out_tree->Write();
		hist->Write();
		out_file->Close();
	}
	
	can1->Close();
	file->Close();	
	rootapp->Delete();
		
	return 0;
}

// For CINT
//int PulseAnalyzer(int argc, char* argv[]){ return main(argc, argv); }
