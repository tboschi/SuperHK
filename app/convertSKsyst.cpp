#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

int main(int argc, char**argv)
{
	TFile *inf  = new TFile(argv[1]);
	TFile *outf = new TFile(argv[2], "RECREATE");

	TTree * syst = (TTree*) inf->Get("sigmatree");

	TString * FijName = new TString;
	double Sigma;

	syst->SetBranchAddress("FijName",  &FijName    );
	syst->SetBranchAddress("Sigma"  ,  &Sigma      );

	int natmo = 2224;

	int index = 0;
	for (int i = 0; i < syst->GetEntries(); ++i)
	{
		syst->GetEntry(i);

		std::string title(FijName->Data());
		if (title.find("t2k") != std::string::npos)
			continue;

		TH1D* h = (TH1D*) inf->Get( FijName->Data() ); 

		std::string name = "systematic_" + std::to_string(index);
		TH1D *newh = new TH1D(name.c_str(), title.c_str(), natmo, 0, natmo);

		for (int n = 1; n < newh->GetNbinsX()+1; ++n)
			newh->SetBinContent(n, h->GetBinContent(n) + 1);

		outf->cd();
		newh->Write();
		++index;
	}

	inf->Close();
	outf->Close();

	return 0;
}
