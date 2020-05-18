#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

#include "tools/CardDealer.h"
#include "physics/ParameterSpace.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[3];
	CardDealer *cd = new CardDealer(cardFile);
	ParameterSpace *parms = new ParameterSpace(cd);

	int kVerbosity;
	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = false;

	//std::map<std::string, std::vector<double> > parms;
	//cd->Get("parm_", parms);

	/*
	for (ib = parms.begin(); ip != parms.end(); ++ip)
	{
		//if less than three elements or fixed parameter (0 or 1 points), skip
		if (ip->second.size() < 3 || ip->second.at(2) < 2)
			continue;
		
		//create array with bins
		std::vector<double> bins(ip->second.at(2) + 1);

		//bin width
		double bw = std::abs(ip->second.at(1) - ip->second.at(0)) / 
			    (ip->second.at(2) - 1);

		//fill first bin
		bins[0] = std::min(ip->second.at(0), ip->second.at(1)) - bw / 2.;

		//and add step until end
		for (int i = 1; i < ip->second.at(2) + 1; ++i)
			bins[i] = bins[i-1] + bw;

		varmap[ip->first] = new double;		//null pointer
		binning[ip->first] = bins;
	}
	*/

	//map of bins for histograms
	ParameterSpace::Binning binning = parms->GetHistogramBinning();
	ParameterSpace::Binning::iterator ix, iy;

	std::map<std::string, double*> varmap;	//map to memory address for variables
	for (ix = binning.begin(); ix != binning.end(); ++ix)
	{
		std::cout << "first " << ix->first << std::endl;
		varmap[ix->first] = new double;		//null pointer
	}


	//map of combination of parameters so that each histogram can be linked
	//to a set of parameters
	std::map<std::string, std::vector<double*> > variables;
	std::map<std::string, TH1*> histograms, histpoints;

	for (ix = binning.begin(); ix != binning.end(); ++ix)
	{
		std::string name = "X2min" + ix->first;
		std::string npnt = "Pmin"  + ix->first;

		std::cout << name << std::endl;

		//create histogram for variable ix->first
		//passing length of binning and address to vector memory
		TH1 *hist = new TH1D(name.c_str(), name.c_str(),
				     binning[ix->first].size()-1, &binning[ix->first][0]);
		TH1 *hpnt = new TH1I(npnt.c_str(), npnt.c_str(),
				     binning[ix->first].size()-1, &binning[ix->first][0]);

		//create vector with pointer to vars
		std::vector<double*> vars;
		vars.push_back(varmap[ix->first]);

		histograms[name] = hist;
		histpoints[name] = hpnt;
		variables[name] = vars;

		//second loop for contouring
		for (iy = std::next(ix); iy != binning.end(); ++iy)
		{
			//these pointers will be used for accessing the map
			ParameterSpace::Binning::iterator tx = ix;
			ParameterSpace::Binning::iterator ty = iy;

			//map keeps the elements sorted, so they are accessed in this order
			//CP, M23, S13, S23, 	and the loop without swapping will be
			//CP vs M23, CP vs S13, CP vs S23, M23 vs S13, M23 vs S23, S13 vs S23
			//but we want this order
			//CP vs M23, S13 vs CP, S23 vs CP, S13 vs M23, S23 vs M23, S13 vs S23
			//so if ix has (C or M) and iy has S swap
			if ( ( (ix->first.find_first_of('C') != std::string::npos) ||
			       (ix->first.find_first_of('M') != std::string::npos) ) &&
			        iy->first.find_first_of('S') != std::string::npos )
			{
				tx = iy;
				ty = ix;
			}

			//create histogram for variable ix->first and iy->first
			//passing length of binning and address to vector memory
			name = "X2" + tx->first + ty->first;
			npnt = "P" + tx->first + ty->first;
			std::cout << name << std::endl;

			TH1 *cont = new TH2D(name.c_str(), name.c_str(),
					      binning[tx->first].size() - 1, &binning[tx->first][0],
					      binning[ty->first].size() - 1, &binning[ty->first][0]);

			TH1 *cpnt = new TH2D(npnt.c_str(), npnt.c_str(),
					      binning[tx->first].size() - 1, &binning[tx->first][0],
					      binning[ty->first].size() - 1, &binning[ty->first][0]);

			vars.clear();
			vars.push_back(varmap[tx->first]);
			vars.push_back(varmap[ty->first]);

			histograms[name] = cont;
			histpoints[name] = cpnt;
			variables[name] = vars;
		}
	}

	double X2;
	int point;

	TFile *inf = new TFile(argv[1], "READ");
	TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));
	sX2->SetBranchAddress("X2",	&X2);
	sX2->SetBranchAddress("Point",	&point);
	sX2->SetBranchAddress("CP",	varmap["CP"]);
	sX2->SetBranchAddress("M12",	varmap["M12"]);
	sX2->SetBranchAddress("M23",	varmap["M23"]);
	sX2->SetBranchAddress("S12",	varmap["S12"]);
	sX2->SetBranchAddress("S13",	varmap["S13"]);
	sX2->SetBranchAddress("S23",	varmap["S23"]);

	//std::string chi2[4] = {"X2minCP", "X2minM23", "X2minS13", "X2minS23"};
	//std::string cont[6] = {"X2CPM23", "X2S13CP", "X2S23CP", "X2S13M23", "X2S23M23", "X2S13S23"};

	std::map<std::string, TH1*>::iterator ih;
	for (int n = 0; n < sX2->GetEntries(); ++n)
	{
		sX2->GetEntry(n);

		for (ih = histograms.begin(); ih != histograms.end(); ++ih)
		{
			int ibin = 0;	//find global bin
			if (variables[ih->first].size() == 1)
				ibin = ih->second->FindBin(*variables[ih->first].at(0));
			else if (variables[ih->first].size() == 2)
				ibin = ih->second->FindBin(*variables[ih->first].at(0),
							   *variables[ih->first].at(1));

			if (ih->second->GetBinContent(ibin) == 0 ||
			    ih->second->GetBinContent(ibin) > X2)
			{
				ih->second->SetBinContent(ibin, X2);
				histpoints[ih->first]->SetBinContent(ibin, point);
			}
		}
	}

	TFile *outf = new TFile(argv[2], "RECREATE");
	outf->cd();

	for (ih = histograms.begin(); ih != histograms.end(); ++ih)
		ih->second->Write();
	for (ih = histpoints.begin(); ih != histpoints.end(); ++ih)
		ih->second->Write();

	outf->Close();
	inf->Close();

	delete cd;

	std::map<std::string, double*>::iterator iv;
	for (iv = varmap.begin(); iv != varmap.end(); ++iv)
		delete iv->second;
	varmap.clear();

	return 0;
}
