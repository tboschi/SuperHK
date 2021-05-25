#include <fstream>
#include <iostream>
#include <chrono>

#include "event/Sample.h"
#include "event/BeamSample.h"
#include "event/ChiSquared.h"

#include "physics/Oscillator.h"
#include "physics/ParameterSpace.h"

#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

//make a likelihood scan - obsx2 as a function of a systematic parameter

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Validation: need one parameters: card" << std::endl;
		return 1;
	}

	// main card
	std::string cardFile = argv[1];
	// main card
	CardDealer cd(cardFile);

	std::string fit_card, osc_card, sample_card;
	if (!cd.Get("oscillation_parameters", osc_card)) {
		std::cerr << "Fitter: no oscillation options card defined, very bad!" << std::endl;
		return 1;
	}
	std::shared_ptr<Oscillator> osc(new Oscillator(osc_card));
	std::unique_ptr<ParameterSpace> parms(new ParameterSpace(osc_card));

	if (!cd.Get("fit_parameters", fit_card)) {
		std::cerr << "Fitter: no fit options card defined, very bad!" << std::endl;;
		return 1;
	}
	std::unique_ptr<ChiSquared> fitter(new ChiSquared(fit_card));


	if (cd.Get("beam_parameters", sample_card))
		fitter->Add<BeamSample>(sample_card, "RBS"); // don't discard empty bins
	if (cd.Get("atmo_parameters", sample_card))
		fitter->Add<AtmoSample>(sample_card);

	// combining samples
	if (!fitter->Combine()) {
		std::cerr << "Fitter: not possible to combine samples. Make sure at least one is defined" << std::endl;
		return 1;
	}

	std::string outfile;
	if (!cd.Get("output", outfile))
		outfile = "validation.dat";
	else if (outfile.find(".dat") == std::string::npos)
		outfile += ".dat";
	
	// parameters for the main script
	int kVerbosity;
	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;

	std::string trueOrder, fitOrder;
	if (!cd.Get("true_hierarchy", trueOrder))
		trueOrder = "normal";
	if (!cd.Get("fit_hierarchy", fitOrder))
		fitOrder = "normal";

	int truePoint, fitPoint;
	if (!cd.Get("true_point", truePoint))
		truePoint = parms->GetNominalEntry();
	if (!cd.Get("fit_point", fitPoint))
		fitPoint = -1;

	Eigen::VectorXd nooscSpectra = fitter->ConstructSamples();

	double tM12, fM12;
	double tM23, fM23;
	double tS12, fS12;
	double tS13, fS13;
	double tS23, fS23;
	double tdCP, fdCP;

	//get True point
	parms->GetEntry(truePoint, tM12, tM23, tS12, tS13, tS23, tdCP);

	if (trueOrder == "normal")
		osc->SetMasses<Oscillator::normal>(tM12, tM23);
	else if (trueOrder == "inverted")
		osc->SetMasses<Oscillator::inverted>(tM12, tM23);
	osc->SetPMNS<Oscillator::sin2>(tS12, tS13, tS23, tdCP);

	//create an object for true spectra (unosc, no systematics), using osc as defined just above
	//i.e. using true_point
	Eigen::VectorXd trueSpectra = fitter->ConstructSamples(osc);

	//create object to store fitSpectra, i.e. using fit_point (no systematics yet)
	Eigen::VectorXd fitSpectra;
	if (fitPoint >= 0) {
		//get Fit point
		parms->GetEntry(fitPoint, fM12, fM23, fS12, fS13, fS23, fdCP);
		if (fitOrder == "normal")
			osc->SetMasses<Oscillator::normal>(fM12, fM23);
		else if (fitOrder == "inverted")
			osc->SetMasses<Oscillator::inverted>(fM12, fM23);
		osc->SetPMNS<Oscillator::sin2>(fS12, fS13, fS23, fdCP);
		fitSpectra = fitter->ConstructSamples(osc);
	}

	//create array to store chi2
	Eigen::ArrayXd chi2;

	//create an array of values of energy scale systematic to generate likelihood scan
	Eigen::VectorXd syst_values = Eigen::VectorXd::LinSpaced(20, -1., 1.);
	std::cout << "systematics vector is " << syst_values << std::endl; 

	//get the number of systematics to define the systematics array
	int num_sytematics = fitter->NumSys();
	//create array for systematics of size num_systematics
	Eigen::VectorXd epsilon = Eigen::VectorXd::Zero(num_sytematics);
	 
	
	//make a copy of true spectra
	Eigen::VectorXd trueSpectra_copy = trueSpectra;
	
	//define output file object
	std::string output = outfile;
	std::ofstream fout(output.c_str());

	//loop over systematic parameter values
	for (int syst_i=0; syst_i<syst_values.size(); ++syst_i){
		//set value for desired systematic
		epsilon(num_sytematics-1) = syst_values(syst_i);

		//store value of systematic parameter
		fout << syst_values(syst_i);

		//vary truespectra, i.e. unoscillated spectra according to systematics set by this epsilon
		chi2 = fitter->ObsX2n(trueSpectra, trueSpectra, epsilon);
		chi2 = chi2.unaryExpr([](double v) { return (std::isfinite(v)||!std::isnan(v))? v : 0.0; });
		std::cout << "chi2 is " << chi2.transpose() << " of size" << chi2.size() <<std::endl;
		Eigen::VectorXd thisSample;

		//loop over samples
		//for now, harcoded for beam samples (4 samples)
		for (int spec_i = 0; spec_i<4; spec_i++){
			 	thisSample = chi2.segment(spec_i*87,87);
				fout << "\t" << thisSample.sum();
				std::cout << "segment is (" << spec_i*87 << "," << (spec_i+1)*87-1 << ")" << std::endl;
				std::cout << "size is " << thisSample.size() <<std::endl;
		}
		fout << "\t" << chi2.sum() << std::endl;
	}
	return 0;
}
