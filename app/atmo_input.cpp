#include <fstream>
#include <iostream>
#include <chrono>

#include "event/AtmoSample.h"

#include "physics/Oscillator.h"
#include "physics/ParameterSpace.h"

#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		std::cerr << "Atmo_input: need three parameters: id cpu card" << std::endl;
		return 1;
	}

	int id  = std::atoi(argv[1]);	//id of current process
	int all = std::atoi(argv[2]);	//total number of processes

	if (id >= all)
	{
		std::cerr << "Atmo_input: requesting process " << id
			  << ", but only " << all << " jobs" << std::endl;
		return 1;
	}

	// main card
	std::string cardFile = argv[3];
	CardDealer cd(cardFile);

	std::string osc_card, sample_card;
	if (!cd.Get("oscillation_parameters", osc_card)) {
		std::cerr << "Atmo_input: no oscillation options card defined, very bad!" << std::endl;
		return 1;
	}
	std::shared_ptr<Oscillator> osc(new Oscillator(cd));
	std::unique_ptr<ParameterSpace> parms(new ParameterSpace(cd));

	if (!cd.Get("atmo_parameters", sample_card)) {
		std::cerr << "Atmo_input: no atmospheric sample card defined, very bad!" << std::endl;
		return 1;
	}
	std::unique_ptr<AtmoSample> as(new AtmoSample(sample_card, "RB"));


	// parameters for the main script
	int kVerbosity;
	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;

	std::string order;
	if (!cd.Get("true_hierarchy", order))
		order = "normal";

	//open output file
	std::string outName;
	cd.Get("output", outName);	
	
	//double *Epsilons = new double[fitter->NumSys()];
	//double *Errors   = new double[fitter->NumSys()];
	double Array[3000];
	int Point, Bins;

	double M12;
	double M23;
	double S12;
	double S13;
	double S23;
	double dCP;

	//get first file for True point

	Eigen::VectorXd spectra = as->ConstructSamples(0);
	std::string atmoArray = "Atmo[" + std::to_string(spectra.size()) + "]/D";

	if (outName.find(".root") == std::string::npos)
		outName += ".root";
	outName.insert(outName.find(".root"), "." + std::to_string(id));

	TFile *outf = new TFile(outName.c_str(), "RECREATE");
	TTree *atmoT = new TTree("atmoTree", "Atmospheric sample");
	atmoT->Branch("Point",	&Point,	"Point/I");
	atmoT->Branch("Bins",	&Bins,	"Bins/I");
	atmoT->Branch("Data",	Array, atmoArray.c_str());

	int entries = parms->GetEntries();

	int off = 0;
	int fpp = entries / all;	//entries per process
	if (id < (entries % all))	//equally distribute
		++fpp;
	else
		off = entries % all;

	int nstart = std::max(off + fpp *  id, 0);
	int nend   = std::min(off + fpp * (id + 1), entries);

	if (kVerbosity)
		std::cout << "Fitter: proc " << id << " / " << all
			  << ", opening files from " << nstart << " to " << nend
			  << " ( " << fpp << " ) " << std::endl;

	for (int i = nstart; i < nend; ++i)
	{
		Point = i;
		parms->GetEntry(Point, M12, M23, S12, S13, S23, dCP);

		if (kVerbosity) {
			std::cout << "\nAtmo_input: now at point " << Point
				  << " (" << i-nstart << "/" << nend-nstart << ")\n";
			std::cout << "m23 " << M23 << ", s13 " << S13
				  << ", s23 " << S23 << ", dcp " << dCP << std::endl;
		}

		if (order == "normal")
			osc->SetMasses<Oscillator::normal>(M12, M23);
		else if (order == "inverted")
			osc->SetMasses<Oscillator::inverted>(M12, M23);
		osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

		spectra = as->ConstructSamples(osc);
		Bins = spectra.size();

		for (int i = 0; i < spectra.size(); ++i)
			Array[i] = spectra(i);

		atmoT->Fill();
	}

	std::cout << "Atmo_input: Finished and out" << std::endl;

	outf->cd();
	atmoT->Write();
	outf->Close();

	return 0;
}
