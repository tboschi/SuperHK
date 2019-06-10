#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

int main(int argc, char** argv)
{
	TFile *f = new TFile(argv[1], "OPEN");
	TFile *o = new TFile(argv[2], "RECREATE");

	TIter next(f->GetListOfKeys());
	TH1D *h = 0;
	TKey *k = 0;
	o->cd();
	while ( k = static_cast<TKey*>(next()) )
	{
		std::string tag(k->GetName());
		if (tag.find_first_of('_') != tag.find_last_of('_'))
			std::cout << "skip " << tag << std::endl;
		else
		{
			h = static_cast<TH1D*>(k->ReadObj());
			h->Write(k->GetName());
		}
	}

	f->Close();
	o->Close();
}
