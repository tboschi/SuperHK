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
		std::cerr << "Vary systematics: need three parameters: card error\n";
		return 1;
	}

	// main card
	std::string cardFile = argv[1];
	// main card
	CardDealer cd(cardFile);

	std::string osc_card, sample_card;
	if (!cd.Get("oscillation_parameters", osc_card)) {
		std::cerr << "Vary systematics: no oscillation options card defined, very bad!\n";
		return 1;
	}
	std::shared_ptr<Oscillator> osc(new Oscillator(osc_card));
	std::unique_ptr<ParameterSpace> parms(new ParameterSpace(osc_card));

	std::unordered_map<std::string, std::shared_ptr<Sample> > samples;
	if (cd.Get("beam_parameters", sample_card))
		samples["beam"] = std::shared_ptr<Sample>(new BeamSample(sample_card, "RBS"));
	if (cd.Get("atmo_parameters", sample_card))
		samples["atmo"] = std::shared_ptr<Sample>(new AtmoSample(sample_card, "RBS"));

	std::string outfile;
	if (!cd.Get("output", outfile))
		outfile = "validation.dat";
	else if (outfile.find(".dat") == std::string::npos)
		outfile += ".dat";
	
	// parameters for the main script
	int kVerbosity;
	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;

	std::string order;
	if (!cd.Get("true_hierarchy", order))
		order = "normal";

	int point;
	if (!cd.Get("true_point", point))
		point = parms->GetNominalEntry();

	double M12, M23, S12, S13, S23, dCP;
	parms->GetEntry(point, M12, M23, S12, S13, S23, dCP);

	if (order == "normal")
		osc->SetMasses<Oscillator::normal>(M12, M23);
	else if (order == "inverted")
		osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	int error = std::strtol(argv[2], NULL, 10);
	//double sigma;
	Eigen::VectorXd sigma = Eigen::VectorXd::LinSpaced(20,-1,1);

	for (const auto &is : samples) {
		if (error >= is.second->NumSys())
			continue;

		// create array for errors
		Eigen::VectorXd epsil = Eigen::VectorXd::Zero(is.second->NumSys());
                // create true spectra using osc as defined just above, i.e. using true_point
                Eigen::VectorXd En = is.second->ConstructSamples(osc);
                auto orig = is.second->Unfold(En);

                for (const auto & ts : orig) {
                        const std::string &type = ts.first;

                        //define output files
                        std::string output = outfile;
                        output.insert(output.find(".dat"), "_" + type);
                        
                        //start filling in output
                        std::ofstream fout(output.c_str());
                        std::cout << "Comparing sample " << type << " and saving to " << output << "\n";
                        
                        //store true values of osc parameters
                        fout << "# true: m12 " << M12 << ", m23 " << M23 << ", s12 " << S12 << ", s13 " << S13
                             << ", s23 " << S23 << ", dcp " << dCP << std::endl;
                        
		for (int syst_i=0;syst_i<sigma.size();++syst_i){
		if (error < 0) {	// want the scale error
			auto scale = is.second->ScaleError();
			for (const auto & se : scale)
				epsil(se.second) = sigma(syst_i);
		}
		else
			epsil(error) = sigma(syst_i); 	// set error in sigma units

		// apply systematics 
		Eigen::VectorXd Ep = is.second->GammaP(En, epsil);

		// poissonian chi2
		Eigen::ArrayXd chi2 = Sample::RawX2n(En, Ep);
		auto x2al = is.second->Unfold(chi2.matrix());
		
		//output
			fout << sigma(syst_i) << "\t" << x2al[type].sum() <<"\n";
		}

		}
	}

	return 0;
}
