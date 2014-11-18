// Optimizer.cpp
// C. Thornsberry
// Nov. 17th, 2014
// Optimize fitting parameters for pulses from a given material
// SYNTAX: ./optimizer {root_filename} {root_treename} {root_branchname} {beta0} {gamma0}

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TApplication.h"

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <time.h>

const double PIXIE_TIME_RES = 4.0; // In ns
const double inv_e = 1.0/std::exp(1.0);

struct linear{
	std::vector<double> x, y;
	
	linear(){}
	
	void clear(){ x.clear(); y.clear(); }
	
	void push_back(double x_, double y_){
		x.push_back(x_); y.push_back(y_);
	}
	
	void Test(){
		unsigned int count = 0;
		std::vector<double>::iterator iterx, itery;
		for(iterx = x.begin(), itery = y.begin(); iterx != x.end() && itery != y.end(); iterx++, itery++){
			std::cout << count << "\t" << *iterx << "\t" << *itery << std::endl;
			count++;
		}
	}
};

static double dabs(const double &input_){
	if(input_ < 0.0){ return -1*input_; }
	return input_;
}

// Mimic the fortran rand() function
double frand(){
	return double(rand())/RAND_MAX;
}

/*static double interpolate(const double &x_, const double &x1_, const double &y1_, const double &x2_, const double &y2_){
	return y2_+(y2_-y1_)*(x_-x1_)/(x2_-x1_);
}*/

static double reverse_interpolate(const double &y_, const double &x1_, const double &y1_, const double &x2_, const double &y2_){
	return -1.0*((y_-y2_)*(x2_-x1_)/(y2_-y1_))+x1_;
}

// 1D function to use for pulse fitting
// x[0] = time t in ns
// parameters: 4
//  par[0] = alpha (normalization)
//  par[1] = phi (phase in ns)
//  par[2] = beta (decay parameter of the exponential in ns)
//  par[3] = gamma (width of the inverted square gaussian in ns^4)
double fit_function(double *x, double *par){
	double arg = x[0] - par[1];
	if(arg >= 0.0){ return par[0]*std::exp(-arg/par[2])*(1 - std::exp(-arg*arg*arg*arg/par[3])); }
	return 0.0;
}

// Calculate the line of best fit (y = alpha + beta*x) for a data set and return the r^2 value
double linear_regression(const std::vector<double> &x_, const std::vector<double> &y_, double &alpha, double &beta){
	double x_mean = 0.0, y_mean = 0.0, xy_mean = 0.0;
	double x2_mean = 0.0, y2_mean = 0.0;
	std::vector<double>::const_iterator iterx, itery;
	
	unsigned int num_values = 0;
	for(iterx = x_.begin(), itery = y_.begin(); iterx != x_.end() && itery != y_.end(); iterx++, itery++){
		x_mean += *iterx;
		y_mean += *itery;
		x2_mean += *iterx*(*iterx);
		y2_mean += *itery*(*itery);
		xy_mean += *iterx*(*itery);
		num_values++;
	}
	
	x_mean = x_mean/num_values;
	y_mean = y_mean/num_values;
	x2_mean = x2_mean/num_values;
	y2_mean = y2_mean/num_values;
	xy_mean = xy_mean/num_values;
	
	double numerator = 0.0, denominator = 0.0;
	for(iterx = x_.begin(), itery = y_.begin(); iterx != x_.end() && itery != y_.end(); iterx++, itery++){
		double temp = (*iterx - x_mean);
		numerator += temp*(*itery - y_mean);
		denominator += temp*temp;
	}
	
	beta = numerator/denominator; 
	alpha = y_mean - beta*x_mean;
	return (xy_mean - x_mean*y_mean)*(xy_mean - x_mean*y_mean)/((x2_mean - x_mean*x_mean)*(y2_mean - y_mean*y_mean));
}

bool process(const std::vector<int> &pulse_, TF1 *function, linear &alpha, linear &phi, double &beta, double &gamma, TCanvas *canvas, bool debug=false){
	if(!function || !canvas){ return false; }

	// Do baseline correction. Use the range to the left of the pulse.
	// The right side of the pulse may have other pulses
	double base_sum = 0.0;
	double darr[pulse_.size()];
	
	unsigned int index = 0;
	unsigned int base_start = pulse_.size() - (unsigned int)(0.05*pulse_.size()); // Take the last 5% of the pulse to do baseline correction
	
	// Convert the pulse integers to double values
	std::vector<int>::const_iterator iter;
	for(iter = pulse_.begin(); iter != pulse_.end(); iter++){
		//darr[index] = *iter + frand(); // Uniformly distribted double value
		darr[index] = (double)(*iter);
		if(index >= base_start){ base_sum += darr[index]; }
		index++;
	}

	// Correct the pulse baseline
	double base_line = base_sum/(pulse_.size()-base_start);
	for(unsigned int i = 0; i < pulse_.size(); i++){
		darr[i] = darr[i] - base_line;
	}

	double *xval = NULL, *yval = NULL;
	TGraph *graph = new TGraph(pulse_.size(), xval, yval);

	// Find the global maximum
	double amplitude = -9999.0;
	unsigned int maximum_bin = 0;
	for(unsigned int i = 0; i < pulse_.size(); i++){
		graph->SetPoint(i, i*PIXIE_TIME_RES, darr[i]);
		if(darr[i] > amplitude){ 
			amplitude = darr[i];
			maximum_bin = i;
		}
	}

	canvas->Clear();
	graph->Draw();
	canvas->Update();
	
	std::string input; std::cout << "\n Use this pulse? (yes/no) "; std::cin >> input;
	if(!(input == "yes" || input == "y")){ return false; }

	// Find the leading edge of the pulse
	unsigned int leading_edge = 0;
	double left, right;
	for(unsigned int i = 1; i < pulse_.size(); i++){
		left = darr[i-1] - 0.1*amplitude; if(left < 0.0){ left = 0.0; }
		right = darr[i] - 0.1*amplitude; if(right < 0.0){ right = 0.0; }
		if((right - left) > 0.0){
			leading_edge = i-1;
			break;
		}
	}
	
	// Find the decay time of the pulse
	double decay_time = 0.0;
	for(unsigned int i = maximum_bin; i < pulse_.size()-1; i++){
		if(darr[i] >= inv_e*amplitude && darr[i+1] <= inv_e*amplitude){
			decay_time = reverse_interpolate(inv_e*amplitude, i*PIXIE_TIME_RES, darr[i], (i+1)*PIXIE_TIME_RES, darr[i+1]) - maximum_bin*PIXIE_TIME_RES;
			break;
		}
	}

	bool locked[4];
	double old_params[4];
	for(unsigned int i = 0; i < 4; i++){ 
		old_params[i] = function->GetParameter(i);
		locked[i] = false;
	}

	double Z = (maximum_bin - leading_edge)*PIXIE_TIME_RES;
	old_params[0] = amplitude*std::exp(Z/decay_time)*(1.0/(1.0-inv_e)); // Alpha (normalization)
	old_params[1] = leading_edge*PIXIE_TIME_RES; // Phi (time offset in ns)
	old_params[2] = decay_time; // Beta (decay time in ns)
	old_params[3] = pow(Z, 4.0); // Gamma (width of square guassian in ns^4);	
	
	function->SetParameter(0, old_params[0]);
	function->SetParameter(1, old_params[1]);
	function->SetParameter(2, old_params[2]);
	function->SetParameter(3, old_params[3]);
	
	if(debug){
		std::cout << " baseline = " << base_line << ", Z = " << Z;
		std::cout << ", amplitude = " << amplitude << ", decay_time = " << decay_time << std::endl; 
		std::cout << "  init[0] = " << old_params[0];
		std::cout << ", init[1] = " << old_params[1];
		std::cout << ", init[2] = " << old_params[2];
		std::cout << ", init[3] = " << old_params[3] << std::endl;
	
		canvas->Clear();
		graph->Draw();
		function->Draw("SAME");
		canvas->Update();
		sleep(1);
	}
	
	unsigned int count = 0;
	bool is_done = false;
	while(!is_done){
		TFitResultPtr func_ptr = graph->Fit(function,"S Q");
		for(unsigned int i = 0; i < 4; i++){ 
			if(locked[i]){ continue; }
			if(dabs(old_params[i] - func_ptr->Value(i)) > 1E-6){ function->SetParameter(i, func_ptr->Value(i)); }
			else{ 
				function->FixParameter(i, func_ptr->Value(i));
				locked[i] = true;
			}
			old_params[i] = function->GetParameter(i);
		}
		
		// Print debug information
		if(debug){ 
			if(locked[0]){ std::cout << "  " << count+1 << ": pars[0] = {" << func_ptr->Value(0) << "}"; }
			else{ std::cout << "  " << count+1 << ": pars[0] = " << func_ptr->Value(0); }
			
			for(unsigned int i = 1; i < 4; i++){ 
				if(locked[i]){ std::cout << ", pars[" << i << "] = {" << func_ptr->Value(i) << "}"; }
				else{ std::cout << ", pars[" << i << "] = " << func_ptr->Value(i); }
			}
			std::cout << ", chi^2 = " << func_ptr->Chi2() << std::endl;
			
			canvas->Clear();
			graph->Draw();
			function->Draw("SAME");
			canvas->Update();
			sleep(1);
		}
		count++;
		
		is_done = true;
		for(unsigned int i = 0; i < 4; i++){ 
			is_done = is_done && locked[i];
			if(!is_done){ break; }
		}
	}

	graph->Delete();

	alpha.push_back(amplitude, old_params[0]);
	phi.push_back(leading_edge*PIXIE_TIME_RES, old_params[1]);
	beta = old_params[2];
	gamma = old_params[3];

	return true;
}

int main(int argc, char* argv[]){
	if(argc < 5){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./optimizer {filename} {treename} {branch} {num_fits} [debug]\n";
		return 1;
	}
	
	bool debug = false;
	if(argc >= 6){
		if(strcmp("debug", argv[5]) == 0){ 
			std::cout << " DEBUGGING...\n";
			debug = true; 
		}
	}
	
	// Seed randomizer
	srand(time(NULL));
	
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");

	// Branch variables
	std::vector<int> wave;
	std::vector<double> energy;
	unsigned int wave_size = 0;
		
	TFile *file = new TFile(argv[1], "READ");
	if(!file->IsOpen()){
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

	TBranch *b_wave;
	tree->SetBranchAddress(argv[3], &wave, &b_wave);
	
	if(!b_wave){
		std::cout << " Failed to load the input branch '" << argv[3] << "'\n";
		file->Close();
		return 1;
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

	TCanvas *can = new TCanvas("can", "canvas");
	can->cd();
	
	double parameters[4];
	parameters[0] = 0.0;
	parameters[1] = 0.0;
	parameters[2] = 0.0;
	parameters[3] = 0.0;
	unsigned int num_fits = (unsigned int)atoi(argv[4]);

	double beta_sum = 0.0;
	double gamma_sum = 0.0;
	linear alpha, phi;
	unsigned int num_fits_done = 0;
	std::cout << " Processing " << tree->GetEntries() << " entries\n";
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		
		TF1 *func = new TF1("func", fit_function, 0.0, wave_size*PIXIE_TIME_RES, 4);
		func->SetLineColor(4); func->SetLineWidth(2);
		func->SetParameters(parameters);	
		
		if(process(wave, func, alpha, phi, parameters[2], parameters[3], can, debug)){ num_fits_done++; }
		beta_sum += parameters[2];
		gamma_sum += parameters[3];
		func->Delete();
		
		if(num_fits_done >= num_fits){ break; }
	}
	
	beta_sum = beta_sum/num_fits;
	gamma_sum = gamma_sum/num_fits;
	
	double A1, B1, A2, B2, R1, R2;
	R1 = linear_regression(alpha.x, alpha.y, A1, B1);
	R2 = linear_regression(phi.x, phi.y, A2, B2);

	std::cout << "\n For t in ns...\n";
	std::cout << "  alpha(t) = " << B1 << "t ";
	if(A1 < 0.0){ std::cout << "- " << -1*A1; }
	else{ std::cout << "+ " << A1; }
	std::cout << ", R^2 = " << R1 << std::endl;

	std::cout << "  phi(t) = " << B2 << "t ";
	if(A2 < 0.0){ std::cout << "- " << -1*A2; }
	else{ std::cout << "+ " << A2; }
	std::cout << ", R^2 = " << R2 << std::endl;

	std::cout << "  beta = " << beta_sum << std::endl;
	std::cout << "  gamma = " << gamma_sum << std::endl;

	can->Close();
	file->Close();	
	rootapp->Delete();
	
	return 0;
}
