#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

class sorter
{
	public:
		sorter(std::vector<double> &v) : vorder(v) {};
		bool operator()(int a, int b)
		{
			return vorder.at(a) < vorder.at(b);
		}
	private:
		std::vector<double> vorder;
};

int main(int argc, char** argv)
{
	//std::map<double, double> compCP, compX2, minCP, minX2;
	std::map<double, double> trueX2, compX2, systX2;

	std::string pen;
	if (argc > 2)
		pen = "_" + std::string(argv[2]);
	std::string cmd = "find " + std::string(argv[1]) +
			  //" -name SpaghettiSens.T2HK_penalised.*.root > .tmp_exclusion";
			  " -name \"SpaghettiSens" + pen + ".T2HK.*.root\" > .tmp_exclusion";
	std::cout << cmd << std::endl;
	system(cmd.c_str());

	std::string file;
	std::ifstream listExclusion(".tmp_exclusion");
	while (std::getline(listExclusion, file))
	//for (int f = 1; f < argc; ++f) 
	{
		std::cout << "opening " << file << std::endl;
		TFile inf(file.c_str(), "READ");
		TTree *t = static_cast<TTree*>(inf.Get("stepX2Tree"));
		double X2, sysX2, dCP, tdCP;
		t->SetBranchAddress("X2",  &X2);
		t->SetBranchAddress("SysX2", &sysX2);
		t->SetBranchAddress("CP",  &dCP);
		t->SetBranchAddress("TCP", &tdCP);

		t->GetEntry(0);

		double min_X2 = X2;
		double min_CP = dCP;

		double comp_X2 = -1;
		double comp_CP = -1;

		std::cout << "looping on " << t->GetEntries() << std::endl;
		for (int i = 0; i < t->GetEntries(); ++i)
		{
			t->GetEntry(i);

			if (std::abs(dCP - tdCP) < 1e-5)
			{
				if (!trueX2.count(tdCP))
					trueX2[tdCP] = X2;
				else if (trueX2[tdCP] > X2)
					trueX2[tdCP] = X2;
			}

			if (1 - std::abs(std::cos(dCP)) < 1e-5)
			{
				if (!compX2.count(tdCP))
				{
					compX2[tdCP] = X2;
					systX2[tdCP] = sysX2;
				}
				else if (compX2[tdCP] > X2)
				{
					compX2[tdCP] = X2;
					systX2[tdCP] = sysX2;
				}
			}
		}

		/*
		if (!compCP.count(tdCP) && comp_CP > 0) 
		{
		       compCP[tdCP] = comp_CP;
		       compX2[tdCP] = comp_X2;
		}
		else if (compX2[tdCP] > comp_X2 && comp_CP > 0)
		{
		       compCP[tdCP] = comp_CP;
		       compX2[tdCP] = comp_X2;
		}

		if (!minCP.count(tdCP)) 
		{
		       minCP[tdCP] = min_CP;
		       minX2[tdCP] = min_X2;
		}
		else if (minX2[tdCP] > min_X2)
		{
		       minCP[tdCP] = min_CP;
		       minX2[tdCP] = min_X2;
		}
		*/

		inf.Close();
	}
	listExclusion.close();

	std::ofstream out("Exclusion.dat");
	std::map<double, double>::iterator im;
	for (im = trueX2.begin(); im != trueX2.end(); ++im)
	{
		out << im->first << "\t"
		    << std::sqrt(std::abs(im->second - compX2[im->first])) << "\t"
		    << im->second << "\t" << systX2[im->first] << "\t"
		    << compX2[im->first] << std::endl;
		   // << "\t" << im->second << "\t" << compX2[im->first] << "\t"
		   // << minCP[im->first] << "\t" << minX2[im->first] << std::endl;
	}

	return 0;
}
