#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

struct syst {
	double eps;
	double err;
};

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cerr << "Name file required." << std::endl;
		return 1;
	}

	std::string cmd = "ls " + std::string(argv[1]) + " > .tmp_epsilon";
	std::cout << cmd << std::endl;
	system(cmd.c_str());

	int point;
	double X2, S13, S23, M23, dCP;
	double tS13, tS23, tM23, tdCP;
	double epsil[120], error[120];

	std::map<double, std::vector<syst> > systdCP;
	std::map<double, std::vector<syst> > systM23;
	std::map<double, std::vector<syst> > systS23;
	std::map<double, std::vector<syst> > systS13;
	std::map<double, std::vector<syst> > systFIT;
	std::map<double, std::vector<syst> >::iterator is;

	std::string file;
	std::ifstream listExclusion(".tmp_epsilon");
	while (std::getline(listExclusion, file))
	{
		std::cout << "opening " << file << std::endl;
		TFile inf(file.c_str(), "READ");
		TTree *t = static_cast<TTree*>(inf.Get("stepX2Tree"));
		t->SetBranchAddress("X2",  &X2);
		t->SetBranchAddress("CP",  &dCP);
		t->SetBranchAddress("S13", &S13);
		t->SetBranchAddress("S23", &S23);
		t->SetBranchAddress("M23", &M23);
		t->SetBranchAddress("TCP",  &tdCP);
		t->SetBranchAddress("TS13", &tS13);
		t->SetBranchAddress("TS23", &tS23);
		t->SetBranchAddress("TM23", &tM23);
		t->SetBranchAddress("Point", &point);
		t->SetBranchAddress("Epsilons", epsil);
		t->SetBranchAddress("Errors", error);

		//std::cout << "looping on " << t->GetEntries() << std::endl;
		for (int i = 0; i < t->GetEntries(); ++i) {
			t->GetEntry(i);

			std::vector<syst> allsyst;
			for (int n = 0; n < 120; ++n) {
				syst s = {epsil[n], error[n]};
				allsyst.push_back(s);
			}

			bool in_dCP = (std::abs(std::abs(dCP) - std::abs(tdCP)) < 1e-9);
			bool in_M23 = (std::abs(M23 - tM23) < 1e-9);
			bool in_S13 = (std::abs(S13 - tS13) < 1e-9);
			bool in_S23 = (std::abs(S23 - tS23) < 1e-9);

			if (in_M23 && in_S13 && in_S23)
				systdCP[dCP] = allsyst;

			if (in_dCP && in_S13 && in_S23)
				systM23[M23] = allsyst;

			if (in_dCP && in_M23 && in_S23)
				systS13[S13] = allsyst;

			if (in_dCP && in_M23 && in_S13)
				systS23[S23] = allsyst;

			if (in_dCP && in_M23 && in_S23 && in_S13)
				systFIT[dCP] = allsyst;
		}

		inf.Close();
	}
	listExclusion.close();

	std::string outname = argv[2];
	std::ofstream out;

	if (outname.find(".dat") == std::string::npos)
		outname += ".dat";

	out.open(outname.c_str());
	for (is = systFIT.begin(); is != systFIT.end(); ++is) {
		for (int n = 0; n < is->second.size(); ++n)
			out << is->first << "\t" << n << "\t"
			    << is->second[n].eps << "\t"
			    << is->second[n].err << std::endl;
		out << std::endl << std::endl;
	}
	out.close();

	std::string thisname = outname;

	thisname.insert(thisname.find(".dat"), "_dCP");
	out.open(thisname.c_str());
	for (is = systdCP.begin(); is != systdCP.end(); ++is) {
		for (int n = 0; n < is->second.size(); ++n)
			out << is->first << "\t" << n << "\t"
			    << is->second[n].eps << "\t"
			    << is->second[n].err << std::endl;
		out << std::endl << std::endl;
	}
	out.close();

	thisname = outname;
	thisname.insert(thisname.find(".dat"), "_M23");
	out.open(thisname.c_str());
	for (is = systM23.begin(); is != systM23.end(); ++is) {
		for (int n = 0; n < is->second.size(); ++n)
			out << is->first << "\t" << n << "\t"
			    << is->second[n].eps << "\t"
			    << is->second[n].err << std::endl;
		out << std::endl << std::endl;
	}
	out.close();

	thisname = outname;
	thisname.insert(thisname.find(".dat"), "_S23");
	out.open(thisname.c_str());
	for (is = systS23.begin(); is != systS23.end(); ++is) {
		for (int n = 0; n < is->second.size(); ++n)
			out << is->first << "\t" << n << "\t"
			    << is->second[n].eps << "\t"
			    << is->second[n].err << std::endl;
		out << std::endl << std::endl;
	}
	out.close();

	thisname = outname;
	thisname.insert(thisname.find(".dat"), "_S13");
	out.open(thisname.c_str());
	for (is = systS13.begin(); is != systS13.end(); ++is) {
		for (int n = 0; n < is->second.size(); ++n)
			out << is->first << "\t" << n << "\t"
			    << is->second[n].eps << "\t"
			    << is->second[n].err << std::endl;
		out << std::endl << std::endl;
	}
	out.close();
	return 0;
}
