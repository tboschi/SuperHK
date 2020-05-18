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
	std::string cardFile = argv[1];
	CardDealer *cd = new CardDealer(cardFile);

	//set up oscillation
	Oscillator *osc = new Oscillator(cd);

	//create chi2 fitter
	ChiSquared *fitter = new ChiSquared(cd);
	fitter->Init();

	//open output file
	std::string outName;
	cd->Get("output", outName);
	
	//preapre TTree reading
	double M12;
	double M23;
	double S12;
	double S13;
	double S23;
	double dCP;

	double epsilold[200];
	double epsilnew[200];
	double epsilerr[200];
	double X2old, X2new;	//min value of X2 at given point
	double Told, Tnew;
	int Point;

	std::string input;
	cd->Get("input", input);

	std::string cmd = "ls " + input + " > .tmp_constrain";
	system(cmd.c_str());
	std::string file;
	std::ifstream listExclusion(".tmp_constrain");

	//get first file for True point
	std::getline(listExclusion, file);
	TFile *inf = new TFile(file.c_str(), "READ");
	if (inf->IsZombie())
	{
		std::cout << "file doesn't exist\n";
		return 0;
	}

	TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));

	sX2->SetBranchAddress("TM12", &M12);
	sX2->SetBranchAddress("TM23", &M23);
	sX2->SetBranchAddress("TS12", &S12);
	sX2->SetBranchAddress("TS13", &S13);
	sX2->SetBranchAddress("TS23", &S23);
	sX2->SetBranchAddress("TCP",  &dCP);

	sX2->GetEntry(0);

	//spectra for true point ( On )
	//it is a concat of all events, E_FHC + E_RHC + M_FHC + M_RHC
	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
	Eigen::VectorXd trueSpectra = fitter->ConstructSpectrum(osc);

	inf->Close();

	TFile *outf = new TFile(outName.c_str(), "RECREATE");
	TTree *cX2 = new TTree("cX2", "my tree");
	cX2->Branch("X2_old", &X2old, "X2old/D");
	cX2->Branch("X2_new", &X2new, "X2new/D");
	cX2->Branch("T_old", &Told, "Told/D");
	cX2->Branch("T_new", &Tnew, "Tnew/D");
	cX2->Branch("E_old", epsilold, "Epsilonsold[200]/D");
	cX2->Branch("E_new", epsilnew, "Epsilonsnew[200]/D");
	cX2->Branch("E_err", epsilerr, "Epsilonserr[200]/D");

	//reset listExclusion steam
	listExclusion.clear();
	listExclusion.seekg(0, std::ios::beg);
	//now loop over all files to get the fit spectra ( E_n )
	while (std::getline(listExclusion, file))
	{
		std::cout << "reading " << file << std::endl;
		inf = new TFile(file.c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));

		sX2->SetBranchAddress("Epsilons", epsilold);
		sX2->SetBranchAddress("X2", &X2old);
		sX2->SetBranchAddress("Time", &Told);
		sX2->SetBranchAddress("Point", &Point);

		sX2->SetBranchAddress("M12", &M12);
		sX2->SetBranchAddress("M23", &M23);
		sX2->SetBranchAddress("S12", &S12);
		sX2->SetBranchAddress("S13", &S13);
		sX2->SetBranchAddress("S23", &S23);
		sX2->SetBranchAddress("CP",  &dCP);

		for (int i = 0; i < sX2->GetEntries(); ++i)
		{
			sX2->GetEntry(i);

			auto t_begin = std::chrono::high_resolution_clock::now();

			//Get expected spectrum ( En )
			osc->SetMasses<Oscillator::normal>(M12, M23);
			osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
			Eigen::VectorXd fitSpectra = fitter->ConstructSpectrum(osc);

			Eigen::VectorXd eps = fitter->FitX2(trueSpectra, fitSpectra);
			Eigen::MatrixXd var = fitter->Covariance(trueSpectra, fitSpectra, eps);

			X2new = fitter->X2(trueSpectra, fitSpectra, eps);

			double diff = 0;
			for (int c = 0; c < eps.size(); ++c)
			{
				epsilnew[67+c] = eps(c);
				epsilerr[67+c] = sqrt(var(c, c));
				diff += pow(epsilold[67+c] - eps(c), 2);
			}
			diff = sqrt(diff);

			std::cout << "\nentry " << i << " computes " << X2new
				  << "  vs  " << X2old << ", dist " << diff << std::endl;
			std::cout << "m23 " << M23 << ", s13 " << S13
				  << ", s23 " << S23 << ", dcp " << dCP << std::endl;
			std::cout << std::endl;


			auto t_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> t_duration = t_end - t_begin;

			Tnew = t_duration.count();

			cX2->Fill();
		}

		inf->Close();
	}

	outf->cd();
	cX2->Write();
	outf->Close();

	delete cd;
	delete osc;
	delete fitter;

	return 0;
}
