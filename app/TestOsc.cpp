#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

#include "physics/Oscillator.h"
#include "probWrapper.h"
#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv)
{
	std::ofstream out;
	TFile *outf;
	if (strcmp(argv[1], "osc") == 0)
	{
		out.open("osctest.dat");
		outf = new TFile("osctest.root", "RECREATE");
	}
	else if (strcmp(argv[1], "prob") == 0)
	{
		out.open("probtest.dat");
		outf = new TFile("probtest.root", "RECREATE");
	}

	double en;
	double oscEE_0, oscMM_0, oscME_0;
	double oscEE_B, oscMM_B, oscME_B;

	TTree *tr = new TTree("osc", "osc");
	tr->Branch("Energy", &en, "en/D");
	tr->Branch("E0_E0", &oscEE_0, "oscEE_0/D");
	tr->Branch("M0_M0", &oscME_0, "oscMM_0/D");
	tr->Branch("M0_E0", &oscME_0, "oscME_0/D");
	tr->Branch("EB_EB", &oscEE_B, "oscEE_B/D");
	tr->Branch("MB_MB", &oscMM_B, "oscMM_B/D");
	tr->Branch("MB_EB", &oscME_B, "oscME_B/D");

	//parameters squared
	double M12 = 7.6e-5;
	double M23 = 2.4e-3;
	double S12 = 0.320;
	double S13 = 0.0256584;
	double S23 = 0.5;
	double dCP = 0;
	//experimental params
	double baseline = 295.0;	//km
	double density = 2.6;		//g/cm³

	std::vector<double> vL(1, baseline), vD(1, density);
	Oscillator *osc = new Oscillator(vL, vD);
	osc->SetMasses<Oscillator::difference>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	ProbWrapper *prob = new ProbWrapper(M23, S23, M12, S13, S12, dCP);
	prob->SetBaseLine(baseline);
	prob->SetDensity(density);

	static const double arr[] = {0.0, 0.05, 0.1, 0.15,
				     0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375,
				     0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575,
				     0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775,
				     0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975,
				     1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45,
				     1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95,
				     2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45,
				     2.5, 2.6, 2.7, 2.8, 2.9,
				     3.0, 3.1, 3.2, 3.3, 3.4,
				     3.5, 3.6, 3.7, 3.8, 3.9,
				     4.0, 4.2, 4.4, 4.6, 4.8,
				     5.0, 5.2, 5.4, 5.6, 5.8,
				     6.0, 6.5, 7.0, 7.5,
				     8, 9, 10, 30};
	std::vector<double> bin(arr, arr + sizeof(arr) / sizeof(arr[0]));

	TH1D* hoscEE_0 = new TH1D("oscEE_0", "oscEE_0", bin.size()-1, arr);
	TH1D* hoscMM_0 = new TH1D("oscMM_0", "oscMM_0", bin.size()-1, arr);
	TH1D* hoscME_0 = new TH1D("oscME_0", "oscME_0", bin.size()-1, arr);
	TH1D* hoscEE_B = new TH1D("oscEE_B", "oscEE_B", bin.size()-1, arr);
	TH1D* hoscMM_B = new TH1D("oscMM_B", "oscMM_B", bin.size()-1, arr);
	TH1D* hoscME_B = new TH1D("oscME_B", "oscME_B", bin.size()-1, arr);

	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < bin.size()-1; ++i)
	{
		en = (bin.at(i+1) + bin.at(i)) / 2.0;
		if (strcmp(argv[1], "osc") == 0)
		{
			oscEE_0 = osc->Probability(Nu::E_, Nu::E_, en);
			oscMM_0 = osc->Probability(Nu::M_, Nu::M_, en);
			oscME_0 = osc->Probability(Nu::M_, Nu::E_, en);
			oscEE_B = osc->Probability(Nu::Eb, Nu::Eb, en);
			oscMM_B = osc->Probability(Nu::Mb, Nu::Mb, en);
			oscME_B = osc->Probability(Nu::Mb, Nu::Eb, en);
		}
		else if (strcmp(argv[1], "prob") == 0)
		{
			oscEE_0 = prob->GetProbNuENuE(en);
			oscMM_0 = prob->GetProbNuMuNuMu(en);
			oscME_0 = prob->GetProbNuMuNuE(en);
			oscEE_B = prob->GetProbNuEBarNuEBar(en);
			oscMM_B = prob->GetProbNuMuBarNuMuBar(en);
			oscME_B = prob->GetProbNuMuBarNuEBar(en);
		}
		tr->Fill();
		hoscEE_0->Fill(en, oscEE_0);
		hoscMM_0->Fill(en, oscMM_0);
		hoscME_0->Fill(en, oscME_0);
		hoscEE_B->Fill(en, oscEE_B);
		hoscMM_B->Fill(en, oscMM_B);
		hoscME_B->Fill(en, oscME_B);
		
		out << bin.at(i) << "\t";
		out << oscEE_0 << "\t" << oscMM_0 << "\t" << oscME_0 << "\t"
		    << oscEE_B << "\t" << oscMM_B << "\t" << oscME_B << std::endl;
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

	if (outf->IsOpen())
	{
		tr->Write();
		hoscEE_0->Write();
		hoscMM_0->Write();
		hoscME_0->Write();
		hoscEE_B->Write();
		hoscMM_B->Write();
		hoscME_B->Write();
		outf->Close();
	}

	return 0;
}
