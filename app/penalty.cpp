#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

void output(std::ofstream &out, std::map<std::string, double> &mInt);
int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	std::string files;
	cd->Get("files", files);

	double M12, sm12 = -1, M23, sm23 = -1;
	double S12, ss12 = -1, S13, ss13 = -1, S23, ss23 = -1;

	std::vector<double> vValue, vSigma;
	std::vector<std::string> vParm;
	std::map<std::string, std::vector<double> > mvals;
	std::map<std::string, std::vector<double> >::iterator im;

	if (cd->Get("p_", mvals))
	{
		for (im = mvals.begin(); im != mvals.end(); ++im)
		{
			vParm.push_back(im->first);
			vValue.push_back(im->second.front());
			vSigma.push_back(im->second.back());
		}
	}

	std::string cmd = "ls " + files + " > .tmp_list";
	system(cmd.c_str());

	std::ifstream in(".tmp_list");
	std::string name;
	while (std::getline(in, name))
	{
		TFile *inf = new TFile(name.c_str(), "OPEN");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		name = name.replace(name.find(".root"), 5, ".dat");
		std::ofstream out(name.c_str());

		TH1D* h = static_cast<TH1D*>(inf->Get("X2minCP"));

		for (int i = 1; i < h->GetNbinsX()+1; ++i)
		{
			double penalty = 0;
			for (int p = 0; p < vParm.size(); ++p)
			{
				std::string hname = "X2minCP_" + vParm[p];
				TH1D *hp = static_cast<TH1D*>(inf->Get(hname.c_str()));

				double chi2 = (hp->GetBinContent(i) - vValue[p]) / vSigma[p];
				std::cout << vParm[p] << "\t" << chi2*chi2 << std::endl;

				penalty += chi2*chi2;
			}

			out << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << "\t" << h->GetBinContent(i) + penalty << "\n";
		}

		inf->Close();
		out.close();
		delete inf;
	}

	in.close();


	return 0;
}
