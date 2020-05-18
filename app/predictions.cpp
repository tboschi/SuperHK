#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);


	double M12, M23;
	double S12, S13, S23;
	double dCP;
	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	std::vector<double> parms;
	parms.push_back(M12);
	parms.push_back(M23);
	parms.push_back(S12);
	parms.push_back(S13);
	parms.push_back(S23);
	parms.push_back(dCP);

	Oscillator *osc = new Oscillator(cd);
	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);

	//reco files	
	//neutrino modes
	std::string reco_file;

	std::string mode[6] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string chan[6] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	std::string horn[2] = {"FHC", "RHC"};
	//std::string name[2] = {texFHC+".tex",  texRHC+".tex"};
	//std::string root[2] = {texFHC+".root", texRHC+".root"};

	Nu fin[6]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu fout[6] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};
	std::string oscf[6] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	Reco* reco;

	std::map<std::string, TH1D*> refSpectra;
	for (int ih = 0; ih < 2; ++ih)	//FHC, RHC
	{
		//rout = new TFile(root[ih].c_str(), "RECREATE");
		for (int im = 0; im < 6; ++im)	//nuE->nuE, nuM->nuE, nuM->nuM
		{
			cd->Get("reco_" + mode[im] + "_" + horn[ih], reco_file);
			reco = new Reco(reco_file);

			const double *bins;
			int nBin = reco->BinsX(bins);
			f1->SetRange(bins[0], bins[nBin-1]);

			std::string fname = horn[ih] + "_" + mode[im];
			TH1D* flux = new TH1D(fname.c_str(), fname.c_str(), nBin, bins);
			//TH1D* flux = static_cast<TH1D*>(oscr->Get(oscf[im].c_str()));
			flux->SetDirectory(0);
			flux->Add(f1);	//filled with ones
			osc->Oscillate(fin[im], fout[im], flux);

			for (int ic = 0; ic < 6; ++ic)	//CCQE, CCnQE, NC x E, M
			{
				//std::cout << "Analysing: " << horn[ih] << ", " << mode[im] << ", " << chan[ic] << std::endl;

				if (chan[ic].find("NC") == std::string::npos)	//it is not a NC
					reco->Scale(chan[ic], flux);
				else if (mode[im] == "nuM0_nuE0" || mode[im] == "nuMB_nuEB") //and no events for these two
					reco->Scale(chan[ic], 0.0);
				//else no need to scale for NC event
					//reco->Scale(chan[ic], flux);

				TH1D* py = reco->Project(chan[ic], 'y');
				if (py)
				{
					std::string hname = chan[ic].substr(0, chan[ic].find_first_of('_'))
							    + "_" + mode[im] + "_" + horn[ih];

					if (refSpectra.count(hname) && refSpectra[hname])
						refSpectra[hname]->Add(py);
					else
						refSpectra[hname] = static_cast<TH1D*>(py->Clone());
				}
			}

			delete flux;
			delete reco;
		}
	}

	int nBins = refSpectra.begin()->second->GetNbinsX();

	std::string base;
	cd->Get("base", base);

	std::ofstream out;

	std::string outfile = base + "/app_FHC.dat";
	out.open(outfile.c_str());
	for (int i = 1; i < nBins+1; ++i)
	{
		out << refSpectra.begin()->second->GetBinLowEdge(i) << "\t"
		    << refSpectra["E_nuM0_nuE0_FHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuMB_nuEB_FHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuE0_nuE0_FHC"]->GetBinContent(i) +
		       refSpectra["E_nuEB_nuEB_FHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuM0_nuM0_FHC"]->GetBinContent(i) +
		       refSpectra["E_nuMB_nuMB_FHC"]->GetBinContent(i) << std::endl;
	}
	out.close();

	outfile = base + "/app_RHC.dat";
	out.open(outfile.c_str());
	for (int i = 1; i < nBins+1; ++i)
	{
		out << refSpectra.begin()->second->GetBinLowEdge(i) << "\t"
		    << refSpectra["E_nuM0_nuE0_RHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuMB_nuEB_RHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuE0_nuE0_RHC"]->GetBinContent(i) +
		       refSpectra["E_nuEB_nuEB_RHC"]->GetBinContent(i) << "\t"
		    << refSpectra["E_nuM0_nuM0_RHC"]->GetBinContent(i) +
		       refSpectra["E_nuMB_nuMB_RHC"]->GetBinContent(i) << std::endl;
	}
	out.close();

	outfile = base + "/disapp_FHC.dat";
	out.open(outfile.c_str());
	for (int i = 1; i < nBins+1; ++i)
	{
		out << refSpectra.begin()->second->GetBinLowEdge(i) << "\t"
		    << refSpectra["M_nuM0_nuM0_FHC"]->GetBinContent(i) << "\t"
		    << refSpectra["M_nuMB_nuMB_FHC"]->GetBinContent(i) << "\t"
		    << refSpectra["M_nuE0_nuE0_FHC"]->GetBinContent(i) +
		       refSpectra["M_nuM0_nuE0_FHC"]->GetBinContent(i) +
		       refSpectra["M_nuEB_nuEB_FHC"]->GetBinContent(i) +
		       refSpectra["M_nuMB_nuEB_FHC"]->GetBinContent(i) << std::endl;
	}
	out.close();

	outfile = base + "/disapp_RHC.dat";
	out.open(outfile.c_str());
	for (int i = 1; i < nBins+1; ++i)
	{
		out << refSpectra.begin()->second->GetBinLowEdge(i) << "\t"
		    << refSpectra["M_nuM0_nuM0_RHC"]->GetBinContent(i) << "\t"
		    << refSpectra["M_nuMB_nuMB_RHC"]->GetBinContent(i) << "\t"
		    << refSpectra["M_nuE0_nuE0_RHC"]->GetBinContent(i) +
		       refSpectra["M_nuM0_nuE0_RHC"]->GetBinContent(i) +
		       refSpectra["M_nuEB_nuEB_RHC"]->GetBinContent(i) +
		       refSpectra["M_nuMB_nuEB_RHC"]->GetBinContent(i) << std::endl;
	}
	out.close();

	delete f1;

	return 0;
}
