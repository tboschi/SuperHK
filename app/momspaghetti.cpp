#include <fstream>
#include <iostream>
#include <chrono>

#include "physics/Oscillator.h"
#include "event/ChiSquared.h"

#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

int main(int argc, char** argv)
{
	int id  = std::atoi(argv[1]);	//id of current process
	int all = std::atoi(argv[2]);	//total number of processes

	std::string cardFile = argv[3];
	CardDealer *cd = new CardDealer(cardFile);

	int kVerbosity;
	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = false;

	std::string trueOrder, fitOrder;
	cd->Get("true_hierarchy", trueOrder);
	cd->Get("fit_hierarchy", fitOrder);

	//set up oscillation
	Oscillator *osc = new Oscillator(cd);

	//create chi2 fitter
	ChiSquared *fitter = new ChiSquared(cd);
	fitter->Init();
	int nsys = fitter->NumSys();

	//open output file
	std::string inName, outName;
	cd->Get("input", inName);	//list of all inputs
	cd->Get("output", outName);	//path for out file with extension
	
	double * Epsilons = new double[fitter->NumSys()];
	double * Errors   = new double[fitter->NumSys()];
	double X2, SysX2;
	double Time;
	int Point, tPoint;

	double M12, tM12;
	double M23, tM23;
	double S12, tS12;
	double S13, tS13;
	double S23, tS23;
	double dCP, tdCP;

	cd->Get("M12", tM12);
	cd->Get("M23", tM23);
	cd->Get("S12", tS12);
	cd->Get("S13", tS13);
	cd->Get("S23", tS23);
	cd->Get("dCP", tdCP);
	cd->Get("Point", tPoint);

	//get first file for True point

	//spectra for true point ( On )
	//it is a concat of all events, E_FHC + E_RHC + M_FHC + M_RHC
	if (trueOrder == "normal")
		osc->SetMasses<Oscillator::normal>(tM12, tM23);
	else if (trueOrder == "inverted")
		osc->SetMasses<Oscillator::inverted>(tM12, tM23);
	osc->SetPMNS<Oscillator::sin2>(tS12, tS13, tS23, tdCP);
	Eigen::VectorXd trueSpectra = fitter->ConstructSpectrum(osc);

	std::string epsilArray = "Epsilons[" +
				 std::to_string(fitter->NumSys()) + "]/D";
	std::string errorArray = "Errors[" +
				 std::to_string(fitter->NumSys()) + "]/D";

	if (outName.find(".root") == std::string::npos)
		outName += ".root";
	outName.insert(outName.find(".root"), "." + std::to_string(id));

	TFile *outf = new TFile(outName.c_str(), "RECREATE");
	TTree *stepX2 = new TTree("stepX2Tree", "Fit Axis Info");
	stepX2->Branch("Time",		&Time, "Time/D");
	stepX2->Branch("Epsilons",	Epsilons, epsilArray.c_str());
	stepX2->Branch("Errors",	Errors,   errorArray.c_str());
	stepX2->Branch("X2",		&X2,		"X2/D");
	stepX2->Branch("SysX2",		&SysX2,		"SysX2/D");
	stepX2->Branch("Point",		&Point,		"Point/D");
	stepX2->Branch("TPoint",	&tPoint,	"TPoint/D");
	stepX2->Branch("CP",	&dCP,	"CP/D");
	stepX2->Branch("TCP",	&tdCP,	"TCP/D");
	stepX2->Branch("M12",	&M12,	"OM12/D");
	stepX2->Branch("TM12",	&tM12,	"TM12/D");
	stepX2->Branch("M23",	&M23,	"M23/D");
	stepX2->Branch("TM23",	&tM23,	"TM23/D");
	stepX2->Branch("S12",	&S12,	"S12/D");
	stepX2->Branch("TS12",	&tS12,	"TS12/D");
	stepX2->Branch("S13",	&S13,	"S13/D");
	stepX2->Branch("TS13",	&tS13,	"TS13/D");
	stepX2->Branch("S23",	&S23,	"S23/D");
	stepX2->Branch("TS23",	&tS23,	"TS23/D");

	std::string file;
	std::ifstream listInput(inName.c_str());	//need to change
	std::vector<std::string> fileInput;

	while (std::getline(listInput, file))
		fileInput.push_back(file);

	int fpp = fileInput.size() / all;		//files per process
	if (id < (fileInput.size() % all))	//equally distribute
		++fpp;

	int nstart = std::max(fpp *  id, 0);
	int nend   = std::min(fpp * (id + 1), int(fileInput.size()));

	if (kVerbosity)
		std::cout << "Proc " << id << " / " << all
			  << ", opening files from " << nstart << " to " << nend << std::endl;

	for (int i = nstart; i < nend; ++i)
	{
		if (kVerbosity)
			std::cout << "reading " << fileInput[i] << std::endl;

		TFile *inf = new TFile(fileInput[i].c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cerr << "file " << fileInput[i] << " doesn't exist\n";
			continue;
		}

		TTree* mc = static_cast<TTree*>(inf->Get("mcTree"));
		if (!mc)
		{
			if (kVerbosity)
				std::cerr << "No mcTree found!" << std::endl;
			continue;
		}


		mc->SetBranchAddress("Point", &Point);
		mc->SetBranchAddress("M12", &M12);
		mc->SetBranchAddress("M23", &M23);
		mc->SetBranchAddress("S12", &S12);
		mc->SetBranchAddress("S13", &S13);
		mc->SetBranchAddress("S23", &S23);
		mc->SetBranchAddress("CP",  &dCP);

		for (int i = 0; i < mc->GetEntries(); ++i)
		{
			mc->GetEntry(i);

			if (kVerbosity)
				std::cout << "\nFitting point " << Point << std::endl;

			auto t_begin = std::chrono::high_resolution_clock::now();

			//Get expected spectrum ( En )
			if (fitOrder == "normal")
				osc->SetMasses<Oscillator::normal>(M12, M23);
			else if (fitOrder == "inverted")
				osc->SetMasses<Oscillator::inverted>(M12, M23);
			osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
			Eigen::VectorXd fitSpectra = fitter->ConstructSpectrum(osc);

			Eigen::VectorXd eps = fitter->FitX2(trueSpectra, fitSpectra);
			Eigen::MatrixXd var = fitter->Covariance(trueSpectra, fitSpectra, eps);

			X2    = fitter->X2(trueSpectra, fitSpectra, eps);
			SysX2 = fitter->SysX2(eps);

			for (int c = 0; c < eps.size(); ++c)	//eps.size == NumSys
			{
				Epsilons[c] = eps(c);
				Errors[c] = sqrt(var(c, c));
			}

			if (kVerbosity)
			{
				std::cout << "m23 " << M23 << ", s13 " << S13
					<< ", s23 " << S23 << ", dcp " << dCP << std::endl;
				std::cout << "X2 computed " << X2 << std::endl;
				std::cout << std::endl;
			}


			auto t_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> t_duration = t_end - t_begin;

			Time = t_duration.count() / 1000.;

			stepX2->Fill();
		}

		inf->Close();
	}

	outf->cd();
	stepX2->Write();
	outf->Close();

	delete cd;
	delete osc;
	delete fitter;

	return 0;
}