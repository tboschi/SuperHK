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
		sorter(std::vector<double> &v) : vorder(v);
	private:
		std::vector<double> vorder;
		bool operator()(int a, int b) return v.at(a) < v.at(b);
};

int main(int argc, char** argv)
{
	std::ifstream inf(argv[1]);
	std::string line;
	std::vector<double> testCP, vcomCP, vcomX2, vminCP, vminX2, exclusion;
	std::vector<int> vI;
	if (false)
	while (std::getline(inf, line))
	{
		TFile f(line.c_str(), "READ");
		std::cout << "opening " << f.GetName() << std::endl;
		TTree *t = static_cast<TTree*>(f.Get("stepX2Tree"));
		double X2, dCP, tdCP;
		t->SetBranchAddress("X2",  &X2);
		t->SetBranchAddress("CP",  &dCP);
		t->SetBranchAddress("TCP", &tdCP);

		t->GetEntry(0);
		testCP.push_back(tdCP);

		double minX2 = X2;
		double minCP = dCP;

		double comX2 = X2;
		double comCP = dCP;

		std::cout << "looping on " << t->GetEntries() << std::endl;
		for (int i = 1; i < t->GetEntries(); ++i)
		{
			t->GetEntry(i);
			if (X2 < minX2)
			{
				minCP = dCP;
				minX2 = X2;
			}

			if (1 - std::abs(std::cos(dCP)) < 1e-5) // it is 0 or ±pi
			{
				if (X2 < comX2)
				{
					comX2 = X2;
					comCP = dCP;
				}
			}
		}

		exclusion.push_back(std::sqrt(std::abs(minX2 - comX2)));
		vcomCP.push_back(comCP);
		vcomX2.push_back(comX2);
		vminCP.push_back(minCP);
		vminX2.push_back(minX2);
		vI.at(i);

		f.Close();
	}

	inf.close();

	std::sort(vI.begin(), vI.end(), sorter(testCP));

	std::ofstream out(argv[2]);
	for (int i = 0; i < testCP.size(); ++i)
	{
		int j = vI.at(i);
		out << testCP.at(j) << "\t" << exclusion.at(j) << "\t"
		    << vcomCP.at(j) << "\t" << vcomX2.at(j) << "\t"
		    << vminCP.at(j) << "\t" << vminX2.at(j) << std::endl;
	}

	return 0;
}
