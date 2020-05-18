#include <iostream>
#include <fstream>

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
		std::cerr << "Need one parameter: card" << std::endl;
		return 1;
	}

	std::string cardFile = argv[1];
	CardDealer *cd = new CardDealer(cardFile);

	int kVerbosity;
	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = false;

	//set up oscillation
	Oscillator *osc = new Oscillator(cd);
	ParameterSpace *parms = new ParameterSpace(cd);

	//create chi2 fitter
	ChiSquared *fitter = new ChiSquared(cd);
	fitter->Init();


	double M12;
	double M23;
	double S12;
	double S13;
	double S23;
	double dCP;


	int truePoint;
	cd->Get("truePoint", truePoint);
	parms->GetEntry(truePoint, M12, M23, S12, S13, S23, dCP);

	std::string trueOrder;
	cd->Get("true_hierarchy", trueOrder);

	if (trueOrder == "normal")
		osc->SetMasses<Oscillator::normal>(M12, M23);
	else if (trueOrder == "inverted")
		osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
	Eigen::VectorXd trueSpectra = fitter->ConstructSpectrum(osc);



	int fitPoint;
	cd->Get("fitPoint", fitPoint);
	parms->GetEntry(fitPoint, M12, M23, S12, S13, S23, dCP);

	std::string fitOrder;
	cd->Get("fit_hierarchy", fitOrder);

	if (fitOrder == "normal")
		osc->SetMasses<Oscillator::normal>(M12, M23);
	else if (fitOrder == "inverted")
		osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
	Eigen::VectorXd fitSpectra  = fitter->ConstructSpectrum(osc);


	if (kVerbosity)
	{
		std::cout << "\nFitting point " << fitPoint << std::endl;
		std::cout << "m23 " << M23 << ", s13 " << S13
			<< ", s24 " << S23 << ", dcp " << dCP << std::endl;
	}

	std::string outName;
	cd->Get("output", outName);	//path for out file with extension
	std::ofstream out(outName.c_str());
	
	for (double s23 = 0.48; s23 < 0.56; s23 += 0.153 / 100)
	//for (double s23 = 0.52; s23 < 0.53; s23 += 0.0001)
	{
		std::cout << "S23 " << s23 << std::endl;
		osc->SetPMNS<Oscillator::sin2>(S12, S13, s23, dCP);
		Eigen::VectorXd fitSpectra  = fitter->ConstructSpectrum(osc);

		Eigen::VectorXd eps = fitter->FitX2(trueSpectra, fitSpectra);
		double SKE = eps(eps.size()-1);
		double bestX2 = fitter->ObsX2(trueSpectra, fitSpectra, eps);

		for (double ske = -1.5; ske < 1.51; ske += 0.01)
		{
			std::cout << "this " << s23 << "\t" << SKE << "\t" << ske << std::endl;
			eps.setZero();
			eps(eps.size()-1) = ske;
			double obsX2 = fitter->ObsX2(trueSpectra, fitSpectra, eps);
			double sysX2 = fitter->SysX2(eps);

			out << s23 << "\t" << SKE << "\t" << bestX2 << "\t" << ske << "\t" << obsX2 << "\t" << sysX2 << std::endl;
		}
		out << "\n" << std::endl;
	}

	delete cd;
	delete osc;
	delete parms;
	delete fitter;

	return 0;
}
