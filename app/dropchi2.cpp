#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "tools/CardDealer.h"

int main(int argc, char** argv)
{
	CardDealer cd(argv[1]);

	std::vector<std::string> allkeys = cd.ListKeys();

	std::vector<std::string> chi2 = {"X2minCP", "X2minM23", "X2minS13", "X2minS23"};
	for (const std::string &x2 : chi2)
	{	
		std::cout << "doing " << x2 << std::endl;

		std::string outAll = x2 + "_all.dat";
		std::string outDif = x2 + "_diff.dat";

		std::ofstream aout(outAll.c_str());
		std::ofstream dout(outDif.c_str());
		
		aout << "#";
		dout << "#";

		std::vector<TH1D*> vh;
		std::vector<TH1I*> vp;
		std::string file;
		for (const std::string k : allkeys) {
			std::string file;
			if (!cd.Get(k, file)) {
				std::cout << "key " << k << " missing file\n";
				continue;
			}
			
			TFile inf(file.c_str(), "READ");
			if (inf.IsZombie()) {
				std::cout << "file " << file << " cannot be opened\n";
				continue;
			}

			if (inf.Get(x2.c_str()))
			{
				TH1D *h = static_cast<TH1D*>(inf.Get(x2.c_str()));
				h->SetDirectory(0);
				vh.push_back(h);
			}

			aout << "\t" << k;
			dout << "\t" << k;
		}
		aout << std::endl;
		dout << std::endl;

		for (int j = 1; j < vh[0]->GetNbinsX()+1; ++j)
		{
			aout << vh[0]->GetBinCenter(j) << "\t" << vh[0]->GetBinContent(j);
			dout << vh[0]->GetBinCenter(j);
			for (size_t i = 1; i < vh.size(); ++i)
			{
				aout << "\t" << vh[i]->GetBinContent(j);
				dout << "\t" << vh[0]->GetBinContent(j) - vh[i]->GetBinContent(j);
			}
			if (vh.size() == 1)
				dout << "\t" << 0;

			aout << std::endl;
			dout << std::endl;
		}

		aout.close();
		dout.close();

		for (size_t i = 0; i < vh.size(); ++i)
			delete vh[i];
	}


	std::vector<std::string> cont = {"X2CPM23", "X2S13CP", "X2S23CP",
					 "X2S13M23", "X2S23M23", "X2S13S23"};
	for (const std::string bd : cont)
	{	
		std::cout << "doing " << bd << std::endl;
		std::string outfile = bd + ".dat";

		std::ofstream out(outfile.c_str());
		
		out << "#x\ty";

		std::vector<TH2D*> vhh;
		for (const std::string k : allkeys) {
			std::string file;
			if (!cd.Get(k, file)) {
				std::cout << "key " << k << " missing file\n";
				continue;
			}
			
			TFile inf(file.c_str(), "READ");
			if (inf.IsZombie()) {
				std::cout << "file " << file << " cannot be opened\n";
				continue;
			}

			if (inf.Get(bd.c_str()))
			{
				TH2D* hh = static_cast<TH2D*>(inf.Get(bd.c_str()));
				hh->SetDirectory(0);
				vhh.push_back(hh);
			}

			out << "\t" << k;
		}
		out << std::endl;

		for (int j = 1; j < vhh[0]->GetXaxis()->GetNbins()+1; ++j)
		{
			for (int k = 1; k < vhh[0]->GetYaxis()->GetNbins()+1; ++k)
			{
				out << vhh[0]->GetXaxis()->GetBinCenter(j) << "\t" << vhh[0]->GetYaxis()->GetBinCenter(k)
					<< "\t" << vhh[0]->GetBinContent(j, k);
				for (size_t i = 1; i < vhh.size(); ++i)
					out << "\t" << vhh[i]->GetBinContent(j, k);
				out << std::endl;
			}
		}

		out.close();

		for (size_t i = 0; i < vhh.size(); ++i)
			delete vhh[i];
	}

	return 0;
}
