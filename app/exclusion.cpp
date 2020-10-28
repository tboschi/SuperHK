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
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " output input [input1 ....]" << std::endl;
		return 1;
	}

	//std::map<double, double> compCP, compX2, minCP, minX2;
	std::map<double, double> trueX2, compX2, systX2;

	//std::map<double, double> minS13, minS23, minM23;
	std::map<double, int> minPoint;
	double S13, S23, M23;
	double TS13, TS23, TM23;
	int point;


	// rest is input files
	for (int f = 2; f < argc; ++f) {
		std::string cmd = "ls " + std::string(argv[f]) + "/*.root > .tmp_exclusion";
		system(cmd.c_str());

		std::string file;
		std::ifstream listExclusion(".tmp_exclusion");
		while (std::getline(listExclusion, file)) {
			TFile inf(file.c_str(), "READ");
			if (inf.IsZombie()) {
				std::cerr << file << " does not exist. Skip\n";
				continue;
			}
			if (!inf.Get("stepX2Tree")) {
				std::cerr << "stepX2Tree not found. Skip\n";
				continue;
			}

			TTree *t = static_cast<TTree*>(inf.Get("stepX2Tree"));
			double X2, sysX2, dCP, tdCP;
			t->SetBranchAddress("X2",  &X2);
			t->SetBranchAddress("SysX2", &sysX2);
			t->SetBranchAddress("CP",  &dCP);
			t->SetBranchAddress("TCP", &tdCP);
			t->SetBranchAddress("S13", &S13);
			t->SetBranchAddress("S23", &S23);
			t->SetBranchAddress("M23", &M23);
			t->SetBranchAddress("TS13", &TS13);
			t->SetBranchAddress("TS23", &TS23);
			t->SetBranchAddress("TM23", &TM23);
			t->SetBranchAddress("Point", &point);

			t->GetEntry(0);

			double min_X2 = X2;
			double min_CP = dCP;

			double comp_X2 = -1;
			double comp_CP = -1;

			if (t->GetEntries() > 0)
				std::cout << "opening " << file
					  << "\twith " << t->GetEntries()
					  << " entries" << std::endl;
			for (int i = 0; i < t->GetEntries(); ++i) {
				t->GetEntry(i);


				// true value of dCP
				if (std::abs(dCP - tdCP) < 1e-5) {
					// never found before or smaller value
					if (!trueX2.count(tdCP) || X2 < trueX2[tdCP])
						trueX2[tdCP] = X2;
				}

				// skip away from test points
				//if (S13 != TS13 || S23 != TS23 || M23 != M23)
					//continue;

				// CP conserving angles
				if (std::abs(std::sin(dCP)) < 1e-5) {
					// never found before or smaller value
					if (!compX2.count(tdCP) || X2 < compX2[tdCP]) {
						compX2[tdCP] = X2;
						//systX2[tdCP] = sysX2;
						minPoint[tdCP] = point;
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
	}

	// first argument is output file
	std::ofstream out(argv[1]);

	for (const auto &im : trueX2) {
		//std::cout << minPoint[im->first] << " -> " << im->first
		//	  << ":\t(" << minS13[im->first]
		//	  << ", " << minS23[im->first] << ", "
		//	  << minM23[im->first] << ")" << std::endl;
		out << im.first << "\t"
		    << std::sqrt(std::abs(im.second - compX2[im.first])) << "\t"
		    << im.second << "\t" //<< systX2[im.first] << "\t"
		    << compX2[im.first] << "\t" << minPoint[im.first] << std::endl;
		   // << "\t" << im->second << "\t" << compX2[im->first] << "\t"
		   // << minCP[im->first] << "\t" << minX2[im->first] << std::endl;
	}

	return 0;
}
