#include <iostream>
#include <iomanip>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
//#include "TIter.h"
#include "TKey.h"
#include "TClass.h"

TH1D* Add(std::map<std::string, TH1D*> &mH, std::string baseName);
int main(int argc, char** argv)
{
	TFile* inF = new TFile(argv[1], "OPEN");
	//TFile* inR = new TFile(argv[2], "OPEN");

	std::map<std::string, TH1D*> mHist;

	TIter next(inF->GetListOfKeys());
	TKey *key;
	std::cout << "Keys found in " << argv[1] << std::endl;
	while ( key = static_cast<TKey*>(next()) )
	{
		if (key->GetClassName() == "TH1D");
		{
			TH1D *h = static_cast<TH1D*>(key->ReadObj());
			h->SetDirectory(0);
			std::cout << "\t> " << h->GetName() << std::endl;
			mHist[std::string(h->GetName())] = h;
		}
	}

	TH1D* hSignal_0_E = Add(mHist, "nuM0_nuE0_E");
	TH1D* hSignal_0_M = Add(mHist, "nuM0_nuE0_M");
	TH1D* hSignal_B_E = Add(mHist, "nuMB_nuEB_E");
	TH1D* hSignal_B_M = Add(mHist, "nuMB_nuEB_M");

	TH1D* hNeutral_E = Add(mHist, "E_NC");
	TH1D* hNeutral_M = Add(mHist, "M_NC");

	TH1D* hElec_0_E = Add(mHist, "nuE0_nuE0_E_CC");
	TH1D* hElec_0_M = Add(mHist, "nuE0_nuE0_M_CC");
	TH1D* hElec_B_E = Add(mHist, "nuEB_nuEB_E_CC");
	TH1D* hElec_B_M = Add(mHist, "nuEB_nuEB_M_CC");

	TH1D* hMuon_0_E = Add(mHist, "nuM0_nuM0_E_CC");
	TH1D* hMuon_0_M = Add(mHist, "nuM0_nuM0_M_CC");
	TH1D* hMuon_B_E = Add(mHist, "nuMB_nuMB_E_CC");
	TH1D* hMuon_B_M = Add(mHist, "nuMB_nuMB_M_CC");


	std::string eRing = std::string(argv[2]) + "_1Re.dat";
	std::string mRing = std::string(argv[2]) + "_1Rm.dat";

	std::ofstream eOut(eRing.c_str());
	std::ofstream mOut(mRing.c_str());

	std::cout << "bins " << mHist.begin()->second->GetNbinsX()+1 << "\n";
	for (int i = 1; i < mHist.begin()->second->GetNbinsX()+1; ++i)
	{
		double en = mHist.begin()->second->GetBinCenter(i);

		eOut << en << "\t";
		mOut << en << "\t";

		eOut << hSignal_0_E->GetBinContent(i) << "\t" << hSignal_B_E->GetBinContent(i) << "\t";
		mOut << hSignal_0_M->GetBinContent(i) << "\t" << hSignal_B_M->GetBinContent(i) << "\t";

		eOut << hNeutral_E->GetBinContent(i) << "\t";
		mOut << hNeutral_M->GetBinContent(i) << "\t";

		eOut << hElec_0_E->GetBinContent(i) << "\t" << hElec_B_E->GetBinContent(i) << "\t";
		mOut << hElec_0_M->GetBinContent(i) << "\t" << hElec_B_M->GetBinContent(i) << "\t";

		eOut << hMuon_0_E->GetBinContent(i) << "\t" << hMuon_B_E->GetBinContent(i) << "\n";
		mOut << hMuon_0_M->GetBinContent(i) << "\t" << hMuon_B_M->GetBinContent(i) << "\n";
	}

	eOut.close();
	mOut.close();

	inF->Close();

	return 0;
}

TH1D* Add(std::map<std::string, TH1D*> &mH, std::string baseName)
{
	TH1D* hTot = 0;
	std::map<std::string, TH1D*>::iterator irh;
	for (irh = mH.begin(); irh != mH.end(); ++irh)
	{
		if (irh->first.find(baseName) != std::string::npos)
		{
			if (!hTot)
			{
				hTot = static_cast<TH1D*>(irh->second->Clone(baseName.c_str()));
				hTot->Reset("ICES");
			}
			hTot->Add(irh->second);
		}
	}
	hTot->SetName(baseName.c_str());
	return hTot;
}
