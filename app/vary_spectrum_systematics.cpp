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

	//create array to store bin content - these are the systematics varied spectra
	Eigen::ArrayXd systVariedSpectra;

	//get the number of systematics to define the systematics array
	int num_sytematics = fitter->NumSys();
	//create array for systematics of size num_systematics
	Eigen::VectorXd epsilon = Eigen::VectorXd::Zero(num_sytematics);
	//set value for desired systematic. here the last entry is energy scale parameter, set to 0.2sigma
	epsilon(num_sytematics-1) = 0.2; 
	
	//make a copy of true spectra
	Eigen::VectorXd trueSpectra_copy = trueSpectra;

	//vary truespectra, i.e. unoscillated spectra accordint to systematics set by epsilon
	systVariedSpectra = fitter->VarySpectrumSystematics(trueSpectra, trueSpectra, epsilon);
	
	//loop over samples
	//for now, harcoded for beam samples (4 samples)
	for (int spec_i = 0; spec_i<4; spec_i++){
		std::string type;
		switch (spec_i) {
			case 0:
				type="E_FHC";
				break;
			case 1:
				type="E_RHC";
				break;
			case 2:
				type="M_RHC";
				break;
			case 3:
				type="M_FHC";
				break;
		}
		//define output files
		std::string output = outfile;
		output.insert(output.find(".dat"), "_" + type);
		//start filling in output
		std::ofstream fout(output.c_str());
		std::cout << "Comparing " << spec_i << " sample and saving to "
			  << output << "\n";

		//store true values of osc parameters
		fout << "# true: m12 " << tM12 << ", m23 " << tM23
		     << ", s12 " << tS12 << ", s13 " << tS13
		     << ", s23 " << tS23 << ", dcp " << tdCP << std::endl;
		if (fitPoint >= 0) {
			//store modified oscillation parameters
			fout << "#  fit: m12 " << fM12 << ", m23 " << fM23
			     << ", s12 " << fS12 << ", s13 " << fS13
			     << ", s23 " << fS23 << ", dcp " << fdCP << std::endl;
		}
		//make a vector to store erec values corresponding to bin number
		std::vector<double> erec_values = fitter->GetErecVect(type);
		//create doubles for storing integrals
		double integral_before = 0;
		double integral_after = 0;
		for (int i = spec_i*87; i < (spec_i+1)*87; ++i) {
			//write energy as first column
			fout << erec_values[i%87] << "\t" << trueSpectra_copy(i);
			//calculate reduced chi2
			double chi2_red = (systVariedSpectra(i)-trueSpectra_copy(i))*(systVariedSpectra(i)-trueSpectra_copy(i))/trueSpectra_copy(i);
			chi2_red = !(std::isnan(chi2_red) || std::isinf(chi2_red)) ? chi2_red : 0;
			//calculate ratio due to systematics
			double ratio = systVariedSpectra(i)/trueSpectra_copy(i);
			ratio = !(std::isnan(ratio) || std::isinf(ratio)) ? ratio : 0;
			if (fitPoint >= 0)
				//store oscillated spectrum (2nd column), systematics varies unosc specrum (3rd column), chi2 (4th column), ratio (5th) column)
				fout << "\t" << fitSpectra(i) <<"\t" << systVariedSpectra(i) << "\t" << chi2_red << "\t" << ratio << "\t" << std::endl;
			else
				fout << "\n";
			//increment integral
			integral_before += trueSpectra_copy(i);
			integral_after += fitSpectra(i);
		}
		//write out integral
		std::cout << "Integral before is " << integral_before << std::endl;
		std::cout << "Integral after is " << integral_after << std::endl;
	
	}

	return 0;
}
