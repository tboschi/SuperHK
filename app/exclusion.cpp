#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "TChain.h"

/*
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
*/

int main(int argc, char** argv)
{
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " output input [input1 ....]" << std::endl;
		return 1;
	}

	//std::map<double, double> compCP, compX2, minCP, minX2;
	std::map<double, double> trueX2, compX2, systX2;
	std::map<double, int> minPoint;

	//std::map<double, double> minS13, minS23, minM23;
	//double S13, S23, M23;
	//double TS13, TS23, TM23;

	// rest is input files
	for (int f = 2; f < argc; ++f) {
		TChain *ch = new TChain("stepX2Tree");
		std::string name(argv[f]);
		name += "/*.root";
		std::cout << "From " << name << "\n";
		std::cout << "\tLoading " << ch->Add(name.c_str()) << " files";

		double X2, dCP, tdCP;
		int point;

		ch->SetBranchAddress("X2",  &X2);
		ch->SetBranchAddress("CP",  &dCP);
		ch->SetBranchAddress("TCP", &tdCP);
		ch->SetBranchAddress("Point", &point);

		ch->GetEntry(0);

		if (ch->GetEntries() > 0)
			std::cout << ": reading " << ch->GetEntries() << " entries";
		for (int i = 0; i < ch->GetEntries(); ++i) {
			ch->GetEntry(i);

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

		std::cout << std::endl;
		delete ch;
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
