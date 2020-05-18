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

void output(std::ofstream &out, double &m12, double &m23, double &s12, double &s13, double &s23, double &dCP);
void output(std::ofstream &out, std::map<std::string, double> &mInt);
int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	TRandom3 *mt = new TRandom3();
	CardDealer *cd = new CardDealer(cardFile);

	std::string texFHC, texRHC;
	cd->Get("out_0", texFHC);
	cd->Get("out_B", texRHC);
	texFHC.erase(texFHC.find(".root"));
	texRHC.erase(texRHC.find(".root"));

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	double random;
	cd->Get("random", random);
	bool krandom = random;

	double M12, M23;
	double S12, S13, S23;
	double dCP;

	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	Oscillator *osc = new Oscillator(densityFile);
	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);
	std::cout << osc->PMNS().cwiseAbs2() << std::endl;

	//std::string oscfile;
	//cd->Get("oscillation", oscfile);
	//TFile *oscr = new TFile(oscfile.c_str());

	//reco files	
	//neutrino modes
	std::string reco_file;

	std::string mode[6] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string chan[6] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	std::string horn[2] = {"FHC", "RHC"};
	std::string name[2] = {texFHC+".tex",  texRHC+".tex"};
	std::string root[2] = {texFHC+".root", texRHC+".root"};

	Nu fin[6]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu fout[6] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};
	std::string oscf[6] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	std::map<std::string, std::string> syserre;
	std::map<std::string, std::string>::iterator is;
	cd->Get("systematic_", syserre);
	std::map<std::string, TH1D*> hsys;
	std::cout << "applying syserre " << syserre.size() << std::endl;
	for (is = syserre.begin(); is != syserre.end(); ++is)
	{
		TFile *sysFile = new TFile(is->second.c_str());
		TIter next(sysFile->GetListOfKeys());
		TKey *k = static_cast<TKey*>(next());

		//each hsys is 1+y, and so summing the first one and then removing 1 from each
		//to have 1 + sum(y)
		std::cout << k->GetName() << std::endl;
		hsys[is->first] = static_cast<TH1D*>(k->ReadObj()->Clone());
		std::cout << hsys[is->first] << std::endl;
		int nBin  = hsys[is->first]->GetNbinsX();
		double lo = hsys[is->first]->GetXaxis()->GetBinLowEdge(1);
		double up = hsys[is->first]->GetXaxis()->GetBinUpEdge(nBin);
		f1->SetRange(lo, up);
		while ( k = static_cast<TKey*>(next()) )
		{
			if (!k->ReadObj()->InheritsFrom("TH1"))
				continue;

			double ran = krandom ? mt->Uniform(-1.0, 1.0) : 1.0;
			TH1D *h1 = static_cast<TH1D*>(k->ReadObj());
			hsys[is->first]->Add(h1,  ran);		//add syst*ran histogram
			hsys[is->first]->Add(f1, -ran);		//remove 1*ran
		}
		hsys[is->first]->SetDirectory(0);
		sysFile->Close();
	}

	std::map<std::string, double> numevts, numsyst;
	Reco* reco;
	TFile *rout;
	for (int ih = 0; ih < 2; ++ih)
	{
		rout = new TFile(root[ih].c_str(), "RECREATE");
		for (int im = 0; im < 6; ++im)
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

			for (int ic = 0; ic < 6; ++ic)
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
					std::string hname = mode[im] + "_" + chan[ic];

					std::string nn;
					if (hname.find("_E_") != std::string::npos)
						nn = "E_" + horn[ih];
					else if (hname.find("_M_") != std::string::npos)
						nn = "M_" + horn[ih];
					
					if (hsys.count(nn))
						py->Multiply(hsys[nn]);

					numevts[hname] = py->Integral();

					rout->cd();
					py->SetName(hname.c_str());
					py->Write(hname.c_str(), TObject::kWriteDelete);
				}
			}

			delete flux;
			delete reco;
		}

		std::ofstream tout(name[ih].c_str());
		tout << "\%" << horn[ih] << " mode\n\n";
		output(tout, M12, M23, S12, S13, S23, dCP);
		output(tout, numevts);
		//output(tout, numsyst);

		tout.close();
		rout->Close();

		numevts.clear();
	}

	delete f1;

	return 0;
}

void output(std::ofstream &out, double &m12, double &m23, double &s12, double &s13, double &s23, double &dCP)
{
	out << std::setfill(' ') << std::setprecision(5) << std::fixed;
	out << "\\textbf{Oscillation parameters}\n\n";
	//1 Ring nu_mu like
	out << "\\begin{tabular}{cccccc}" << "\n";
	out << "\t\\toprule" << "\n";
	out << "\t& $\\Delta m^2{12} / {e-5}\\text{eV}^2$"
	    << "\t& $\\Delta m^2{23} / {e-3}\\text{eV}^2$"
	    << "\t& $\\sin^2 2\\theta_{12}$"
	    << "\t& $\\sin^2 2\\theta_{13}$"
	    << "\t& $\\sin^2  \\theta_{23}$"
	    << "\t& $\\delta_\\text{CP}$"
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";
	//CCQE numbers
	out << std::setw(4) << m12 / 1e-5 << "\t&"
	    << std::setw(4) << m23 / 1e-4 << "\t&"
	    << std::setw(4) << 0.5*(1-sqrt(s12)) << "\t&"
	    << std::setw(4) << 0.5*(1-sqrt(s13)) << "\t&"
	    << std::setw(4) << s23 << "\t&"
	    << std::setw(4) << dCP
	    << "\t\\\\" << "\n";
	out << "\t\\bottomrule" << "\n";
	out << "\\end{tabular}" << "\n\n\n";
}

void output(std::ofstream &out, std::map<std::string, double> &mInt)
{
	out << std::setfill(' ') << std::setprecision(5) << std::fixed;
	out << "\\textbf{1 Ring $\\nu_\\mu$-like}\n\n";
	//1 Ring nu_mu like
	out << "\\begin{tabular}{lrrrrrrr}" << "\n";
	out << "\t\\toprule" << "\n";
	out << "\t& $\\nu_\\mu$"
	     << "\t& $\\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu$"
	     << "\t& $\\cj{\\nu}_e$"
	     << "\t& $\\nu_\\mu \\to \\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu \\to \\cj{\\nu}_e$"
	     << "\t& Total"
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";
	//CCQE numbers
	out << "\tCCQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_M_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCQE"] + mInt["nuE0_nuE0_M_CCQE"] +
				mInt["nuMB_nuMB_M_CCQE"] + mInt["nuEB_nuEB_M_CCQE"] +
				mInt["nuM0_nuE0_M_CCQE"] + mInt["nuMB_nuEB_M_CCQE"]
	     << "\t\\\\" << "\n";
	//CCnQE numbers
	out << "\tCCnQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_M_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCnQE"] + mInt["nuE0_nuE0_M_CCnQE"] +
				mInt["nuMB_nuMB_M_CCnQE"] + mInt["nuEB_nuEB_M_CCnQE"] +
				mInt["nuM0_nuE0_M_CCnQE"] + mInt["nuMB_nuEB_M_CCnQE"]
	     << "\t\\\\" << "\n";
	//NC numbers
	out << "\tNC\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_M_NC"] + mInt["nuE0_nuE0_M_NC"] +
				mInt["nuMB_nuMB_M_NC"] + mInt["nuEB_nuEB_M_NC"] +
				mInt["nuM0_nuE0_M_NC"] + mInt["nuMB_nuEB_M_NC"]
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";	
	//Total numbers
	out << "\tTot\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCQE"] + mInt["nuM0_nuM0_M_CCnQE"] + mInt["nuM0_nuM0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_M_CCQE"] + mInt["nuE0_nuE0_M_CCnQE"] + mInt["nuE0_nuE0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_M_CCQE"] + mInt["nuMB_nuMB_M_CCnQE"] + mInt["nuMB_nuMB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_M_CCQE"] + mInt["nuEB_nuEB_M_CCnQE"] + mInt["nuEB_nuEB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_M_CCQE"] + mInt["nuM0_nuE0_M_CCnQE"] + mInt["nuM0_nuE0_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_M_CCQE"] + mInt["nuMB_nuEB_M_CCnQE"] + mInt["nuMB_nuEB_M_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_M_CCQE"] + mInt["nuM0_nuM0_M_CCnQE"] + mInt["nuM0_nuM0_M_NC"] +
	     			mInt["nuE0_nuE0_M_CCQE"] + mInt["nuE0_nuE0_M_CCnQE"] + mInt["nuE0_nuE0_M_NC"] +
	     			mInt["nuMB_nuMB_M_CCQE"] + mInt["nuMB_nuMB_M_CCnQE"] + mInt["nuMB_nuMB_M_NC"] +
	     			mInt["nuEB_nuEB_M_CCQE"] + mInt["nuEB_nuEB_M_CCnQE"] + mInt["nuEB_nuEB_M_NC"] +
	     			mInt["nuM0_nuE0_M_CCQE"] + mInt["nuM0_nuE0_M_CCnQE"] + mInt["nuM0_nuE0_M_NC"] +
	     			mInt["nuMB_nuEB_M_CCQE"] + mInt["nuMB_nuEB_M_CCnQE"] + mInt["nuMB_nuEB_M_NC"] 
	     << "\t\\\\" << "\n";
	out << "\t\\bottomrule" << "\n";	
	out << "\\end{tabular}" << "\n";

	out << "\n\\vfill\n";

	out << "\\textbf{1 Ring $\\nu_e$-like}\n\n";
	//1 Ring nu_e like
	out << "\\begin{tabular}{lrrrrrrr}" << "\n";
	out << "\t\\toprule" << "\n";
	out << "\t& $\\nu_\\mu$"
	     << "\t& $\\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu$"
	     << "\t& $\\cj{\\nu}_e$"
	     << "\t& $\\nu_\\mu \\to \\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu \\to \\cj{\\nu}_e$"
	     << "\t& Total"
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";
	//CCQE numbers
	out << "\tCCQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_E_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCQE"] + mInt["nuE0_nuE0_E_CCQE"] +
				mInt["nuMB_nuMB_E_CCQE"] + mInt["nuEB_nuEB_E_CCQE"] +
				mInt["nuM0_nuE0_E_CCQE"] + mInt["nuMB_nuEB_E_CCQE"]
	     << "\t\\\\" << "\n";
	//CCnQE numbers
	out << "\tCCnQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_E_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCnQE"] + mInt["nuE0_nuE0_E_CCnQE"] +
				mInt["nuMB_nuMB_E_CCnQE"] + mInt["nuEB_nuEB_E_CCnQE"] +
				mInt["nuM0_nuE0_E_CCnQE"] + mInt["nuMB_nuEB_E_CCnQE"]
	     << "\t\\\\" << "\n";
	//NC numbers
	out << "\tNC\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_E_NC"] + mInt["nuE0_nuE0_E_NC"] +
				mInt["nuMB_nuMB_E_NC"] + mInt["nuEB_nuEB_E_NC"] +
				mInt["nuM0_nuE0_E_NC"] + mInt["nuMB_nuEB_E_NC"]
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";	
	//Total numbers
	out << "\tTot\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCQE"] + mInt["nuM0_nuM0_E_CCnQE"] + mInt["nuM0_nuM0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_E_CCQE"] + mInt["nuE0_nuE0_E_CCnQE"] + mInt["nuE0_nuE0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_E_CCQE"] + mInt["nuMB_nuMB_E_CCnQE"] + mInt["nuMB_nuMB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_E_CCQE"] + mInt["nuEB_nuEB_E_CCnQE"] + mInt["nuEB_nuEB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_E_CCQE"] + mInt["nuM0_nuE0_E_CCnQE"] + mInt["nuM0_nuE0_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_E_CCQE"] + mInt["nuMB_nuEB_E_CCnQE"] + mInt["nuMB_nuEB_E_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_E_CCQE"] + mInt["nuM0_nuM0_E_CCnQE"] + mInt["nuM0_nuM0_E_NC"] +
	     			mInt["nuE0_nuE0_E_CCQE"] + mInt["nuE0_nuE0_E_CCnQE"] + mInt["nuE0_nuE0_E_NC"] +
	     			mInt["nuMB_nuMB_E_CCQE"] + mInt["nuMB_nuMB_E_CCnQE"] + mInt["nuMB_nuMB_E_NC"] +
	     			mInt["nuEB_nuEB_E_CCQE"] + mInt["nuEB_nuEB_E_CCnQE"] + mInt["nuEB_nuEB_E_NC"] +
	     			mInt["nuM0_nuE0_E_CCQE"] + mInt["nuM0_nuE0_E_CCnQE"] + mInt["nuM0_nuE0_E_NC"] +
	     			mInt["nuMB_nuEB_E_CCQE"] + mInt["nuMB_nuEB_E_CCnQE"] + mInt["nuMB_nuEB_E_NC"] 
	     << "\t\\\\" << "\n";
	out << "\t\\bottomrule" << "\n";	
	out << "\\end{tabular}" << "\n\n\n";
}
