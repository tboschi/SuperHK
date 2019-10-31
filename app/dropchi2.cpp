#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TH1.h"
#include "TFile.h"

int main(int argc, char** argv)
{
	std::cout << "Reading " << argc-1 << " files" << std::endl;

	TH1D *h0 = 0;
	std::string chi2[4] = {"X2minCP", "X2minM23", "X2minS13", "X2minS23"};
	for (int x = 0; x < 4; ++x)
	{	
		std::string outAll = chi2[x] + "_all.dat";
		std::string outDif = chi2[x] + "_diff.dat";

		std::ofstream aout(outAll.c_str());
		std::ofstream dout(outDif.c_str());
		
		aout << "#";
		dout << "#";

		std::vector<TH1D*> vh;
		for (int f = 1; f < argc; ++f)
		{
			TFile inf(argv[f], "OPEN");
			if (inf.IsZombie())
				continue;
			TH1D* h = static_cast<TH1D*>(inf.Get(chi2[x].c_str()));
			h->SetDirectory(0);

			vh.push_back(h);
		
			std::string name(argv[f]);
			name.erase(0, name.find("errorstudy/")+11);
			name.erase(name.find_first_of('/'));
			aout << "\t" << name;
			dout << "\t" << name;
		}
		aout << std::endl;
		dout << std::endl;

		for (int j = 1; j < vh[0]->GetNbinsX()+1; ++j)
		{
			aout << vh[0]->GetBinCenter(j) << "\t" << vh[0]->GetBinContent(j);
			dout << vh[0]->GetBinCenter(j);
			for (int i = 1; i < vh.size(); ++i)
			{
				aout << "\t" << vh[i]->GetBinContent(j);
				dout << "\t" << vh[0]->GetBinContent(j) - vh[i]->GetBinContent(j);
			}
			aout << std::endl;
			dout << std::endl;
		}

		aout.close();
		dout.close();

		for (int i = 0; i < vh.size(); ++i)
			delete vh[i];
	}

	return 0;
}
