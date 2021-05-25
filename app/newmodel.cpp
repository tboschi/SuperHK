#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TKey.h"
#include "TList.h"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Addmatrix: need card file\n";
		return 1;
	}

	CardDealer cd(argv[1]);

	std::map<std::string, std::string> matx_files;
	std::map<std::string, std::vector<std::string> > syst_files;

	if (!cd.Get("matx_", matx_files)) {
		std::cerr << "expecting a \"matrices\" key in card file\n";
		return 1;
	}

	if (!cd.Get("syst_", syst_files)) {
		std::cerr << "expecting a \"matrices\" key in card file\n";
		return 1;
	}

	std::set<int> skip_sys;
	cd.Get("skip", skip_sys);	// errors to skip

	std::map<std::string, TMatrixT<double>*> matrices;
	for (const auto &ff : matx_files) {
		TFile inf(ff.second.c_str(), "READ");
		if (inf.IsZombie())
			continue;
		std::cout << "opening " << ff.first << " matrix\n";

		matrices[ff.first] = static_cast<TMatrixT<double>*> (inf.Get("covariance"));
	}

	std::cout << "Collected " << matrices.size() << " matrices" << std::endl;

	if (!syst_files.count("old")) {
		std::cerr << "expecting old systematic files\n";
		return 1;
	}

	if (!syst_files.count("new")) {
		std::cerr << "expecting new systematic files\n";
		return 1;
	}

	std::cout << "Number of files are " << syst_files["old"].size() << std::endl;
	for (int i = 0; i < syst_files["old"].size(); ++i) {
		// it is sample type name

		TFile *sysf = new TFile(syst_files["old"][i].c_str(), "READ");
		if (sysf->IsZombie()) {
			std::cerr << "WARNING - file " << syst_files["old"].at(i)
				<< " does not exist" << std::endl;
			continue;
		}

		TFile *newf = new TFile(syst_files["new"][i].c_str(), "RECREATE");
		if (newf->IsZombie()) {
			std::cerr << "WARNING - file " << syst_files["new"][i]
				<< " does not exist" << std::endl;
			continue;
		}

		TIter next(sysf->GetListOfKeys());
		TH1D* hsys;
		TKey *k;
		int off = 0, k_pre = -1;
		while (k = static_cast<TKey*>(next()))
		{
			std::string sysname = k->GetName();
			std::string title = k->GetTitle();
			if (title.find("= ") != std::string::npos)
				title.erase(title.find("=")+2);
			hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

			std::string postfix;
			if (sysname.find_first_of('_') != sysname.find_last_of('_')) {
				postfix = sysname.substr(sysname.find_last_of('_'));
				sysname.erase(sysname.find_last_of('_'));
				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number now
			}
			else
				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is already a number

			int k_err = std::stoi(sysname);
			if (skip_sys.find(k_err) != skip_sys.end()) {
				if (k_pre != k_err)
					++off;
				k_pre = k_err;
				continue;
			}

			double olderr = sqrt(matrices["old"]->operator()(k_err, k_err));
			double newerr = sqrt(matrices["new"]->operator()(k_err-off, k_err-off));


			TH1D *hnew = static_cast<TH1D*>(hsys->Clone());
			hnew->Reset("ICES");
			for (int n = 1; n < hsys->GetNbinsX()+1; ++n) {
				double cont = (hsys->GetBinContent(n) - 1) * newerr / olderr + 1;
				hnew->SetBinContent(n, cont);
			}

			newf->cd();

			sysname = "syserre_" + std::to_string(k_err-off) + postfix;

			std::cout << "scaling " << k->GetName() << " by " << newerr / olderr
			          << " and saving as " << sysname << std::endl;

			std::stringstream ttl;
			ttl << title << std::setprecision(3) << newerr;
			hnew->SetTitle(ttl.str().c_str());
			hnew->Write(sysname.c_str());
		}

		sysf->Close();
		newf->Close();
	}

	return 0;
}
