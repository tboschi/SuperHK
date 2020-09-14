#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

int main(int argc, char** argv)
{
	std::cout << "Reading " << argc-1 << " files" << std::endl;

	TH1D *h0 = 0;
	std::string chi2[4] = {"X2minCP", "X2minM23", "X2minS13", "X2minS23"};
	std::string pnt2[4] = {"X2minCP", "X2minM23", "X2minS13", "X2minS23"};
	for (int x = 0; x < 4; ++x)
	{	
		std::cout << "doing " << chi2[x] << std::endl;

		std::string outAll = chi2[x] + "_all.dat";
		std::string outDif = chi2[x] + "_diff.dat";

		std::ofstream aout(outAll.c_str());
		std::ofstream dout(outDif.c_str());
		
		aout << "#";
		dout << "#";

		std::vector<TH1D*> vh;
		std::vector<TH1I*> vp;
		for (int f = 1; f < argc; ++f)
		{
			TFile inf(argv[f], "OPEN");
			if (inf.IsZombie())
				continue;

			if (inf.Get(chi2[x].c_str()))
			{
				TH1D *h = static_cast<TH1D*>(inf.Get(chi2[x].c_str()));
				h->SetDirectory(0);
				vh.push_back(h);
			}
			//if (inf.Get(pnt2[x].c_str()))
			//{
			//	TH1I *p = static_cast<TH1I*>(inf.Get(chi2[x].c_str()));
			//	p->SetDirectory(0);
			//	vp.push_back(p);
			//}

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
			if (vh.size() == 1)
				dout << "\t" << 0;

			aout << std::endl;
			dout << std::endl;
		}

		aout.close();
		dout.close();

		for (int i = 0; i < vh.size(); ++i)
			delete vh[i];
	}


	//std::string cont[6] = {"X2CPM23", "X2CPS13", "X2CPS23", "X2M23S13", "X2M23S23", "X2S13S23"};
	std::string cont[6] = {"X2CPM23", "X2S13CP", "X2S23CP", "X2S13M23", "X2S23M23", "X2S13S23"};
	for (int x = 0; x < 6; ++x)
	{	
		std::cout << "doing " << cont[x] << std::endl;
		std::string outfile = cont[x] + ".dat";

		std::ofstream out(outfile.c_str());
		
		out << "#x\ty";

		std::vector<TH2D*> vhh;
		for (int f = 1; f < argc; ++f)
		{
			TFile inf(argv[f], "OPEN");
			if (inf.IsZombie())
				continue;
			TH2D* hh = static_cast<TH2D*>(inf.Get(cont[x].c_str()));
			hh->SetDirectory(0);

			vhh.push_back(hh);
		
			std::string name(argv[f]);
			name.erase(0, name.find("errorstudy/")+11);
			name.erase(name.find_first_of('/'));
			out << "\t" << name;
		}
		out << std::endl;

		for (int j = 1; j < vhh[0]->GetXaxis()->GetNbins()+1; ++j)
		{
			for (int k = 1; k < vhh[0]->GetYaxis()->GetNbins()+1; ++k)
			{
				out << vhh[0]->GetXaxis()->GetBinCenter(j) << "\t" << vhh[0]->GetYaxis()->GetBinCenter(k)
					<< "\t" << vhh[0]->GetBinContent(j, k);
				for (int i = 1; i < vhh.size(); ++i)
					out << "\t" << vhh[i]->GetBinContent(j, k);
				out << std::endl;
			}
		}

		out.close();

		for (int i = 0; i < vhh.size(); ++i)
			delete vhh[i];
	}

	return 0;
}
