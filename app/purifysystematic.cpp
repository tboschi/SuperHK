#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TKey.h"

int main(int argc, char** argv)
{
	TFile *f = new TFile(argv[1], "OPEN");
	TFile *o = new TFile(argv[2], "RECREATE");

	std::vector<std::string> saveList;
	std::vector<std::string>::reverse_iterator is;
	TIter next(f->GetListOfKeys());
	TH1D *h = 0;
	TKey *k = 0;
	o->cd();
	while ( k = static_cast<TKey*>(next()) )
	{
		std::string tag(k->GetName());
		//std::cout << "tag " << tag << std::endl;

		if (tag.find("_p1") != std::string::npos) //remove full systematic and save spline
		{
			//std::cout << "looking for " << tag.substr(0, tag.find("_p1")) << std::endl;
			for (is = saveList.rbegin(); is != saveList.rend(); ++is)
			{
				//std::cout << "\t" << *is << std::endl;
				if (*is == tag.substr(0, tag.find("_p1")))
				{
					//std::cout << "erasing " << *(std::next(is).base()) << std::endl;
					saveList.erase(std::next(is).base());
					break;
				}
			}
			//std::cout << "skip " << *is << std::endl;
		}

		saveList.push_back(tag);
	}

	//std::cout << "Found " << saveList.size() << " systematics" << std::endl;

	for (int i = 0; i < saveList.size(); ++i)
	{
		//std::cout << "saveing " << saveList[i] << std::endl;
		h = static_cast<TH1D*>(f->Get(saveList[i].c_str()));
		h->Write(saveList[i].c_str());
	}

	f->Close();
	o->Close();
}
