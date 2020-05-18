#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "physics/ParameterSpace.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[2];
	CardDealer *cd = new CardDealer(cardFile);
	ParameterSpace *parms = new ParameterSpace(cd);
	ParameterSpace::Binning binning = parms->GetBinning();
	ParameterSpace::Binning::iterator ib;

	std::string phist[4] = {"PminCP", "PminM23", "PminS13", "PminS23"};
	for (int x = 0; x < 4; ++x)
	{	
		std::cout << "reading " << phist[x] << std::endl;

		std::string outAll = phist[x] + "_info.dat";

		std::ofstream aout(outAll.c_str());
		
		TFile inf(argv[1], "OPEN");
		if (inf.IsZombie())
			continue;

		TH1I *p = 0;
		if (inf.Get(phist[x].c_str()))
			p = static_cast<TH1I*>(inf.Get(phist[x].c_str()));
		else
			continue;

		p->SetDirectory(0);

		std::string name(argv[1]);
		name.erase(0, name.find("errorstudy/")+11);
		name.erase(name.find_first_of('/'));
		aout << "#" << name << "\tpoint";
		for (ib = binning.begin(); ib != binning.end(); ++ib)
			aout << "\t" << ib->first;
		aout << std::endl;

		for (int j = 1; j < p->GetNbinsX()+1; ++j)
		{
			aout << p->GetBinCenter(j) << "\t" << p->GetBinContent(j);

			std::map<std::string, double> best = parms->GetEntry(p->GetBinContent(j));
			std::map<std::string, double>::iterator it;
			for (it = best.begin(); it != best.end(); ++it)
				aout << "\t" << it->second;
			aout << std::endl;
		}

		aout.close();
	}

	return 0;
}
