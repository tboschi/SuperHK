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

void beginDocument(std::ofstream &out);
void endDocument(std::ofstream &out);
void printInfo(std::ofstream &out, std::string title, double &m12, double &m23, double &s12, double &s13, double &s23, double &dCP);
void addTable(std::ofstream &out, std::string title, std::map<std::string, double> &mInt);
void addTable(std::ofstream &out, std::string title, std::string p, std::map<std::string, double> &mInt);
int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	TRandom3 *mt = new TRandom3();
	CardDealer *cd = new CardDealer(cardFile);

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	double M12, M23;
	double S12, S13, S23;
	double dCP;
	std::string baseName;

	cd->Get("name", baseName);
	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	Oscillator *osc = new Oscillator(densityFile);
	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);
	//std::cout << osc->PMNS().cwiseAbs2() << std::endl;

	std::string baseOut;
	cd->Get("base", baseOut);
	std::string texFile = baseOut + baseName + "_number_events.tex";
	std::ofstream tout(texFile.c_str());

	beginDocument(tout);

	printInfo(tout, baseName, M12, M23, S12, S13, S23, dCP);


	//reco files	
	//neutrino modes
	std::string reco_file;

	std::string mode[6] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string chan[6] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	std::string horn[2] = {"FHC", "RHC"};

	Nu fin[6]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu fout[6] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};
	std::string oscf[6] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	std::map<std::string, double> numevts, numsyst;
	Reco* reco;
	//TFile *rout;

	//TH1D* FCH1Re, *FHC1Rmu, *RHC1Re, *RHC1Rmu;
	//std::map<std::string, TH1D*> refSpectra;
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
					std::string hname = mode[im] + "_" + chan[ic];
					std::cout << "calculating events for " << hname << std::endl;

					numevts[hname] = py->Integral();
				}
			}

			delete flux;
			delete reco;
		}

		addTable(tout, horn[ih], numevts);

		numevts.clear();
	}

	endDocument(tout);
	tout.close();

	chdir(baseOut.c_str());
	std::string cmd = "pdflatex " + baseName + "_number_events.tex";
	system(cmd.c_str());

	delete f1;

	return 0;
}

void beginDocument(std::ofstream &out)
{
	out << "\\documentclass[aspectratio=169]{beamer}\n";
	out << "\\usepackage{booktabs}\n\n";
	out << "\\setbeamertemplate{navigation symbols}{}\n";
	out << "\\setbeamerfont{frametitle}{size=\\small}\n";
	out <<"\n\\begin{document}\n\n";
}

void endDocument(std::ofstream &out)
{
	out <<"\n\\end{document}\n\n";
}

void printInfo(std::ofstream &out, std::string title, double &m12, double &m23, double &s12, double &s13, double &s23, double &dCP)
{
	out << "\\begin{frame}\n";
	out << "\t\\frametitle{Oscillation parameters}\n";
	out << "\t\\footnotesize\n";
	out << "\t\\centering\n";

	out << std::setfill(' ') << std::setprecision(5) << std::fixed;
	out << "\\textbf{" << title << "}\n\n\\bigskip";

	out << "\\begin{tabular}{cc}\n";
	out << "\t\\toprule\n";
	out << "\tParameter\t&Value\t\\\\\n";
	out << "\t\\midrule\n";
	out << "\t$\\Delta m_{12}^2/ 10^{-5}\\,\\text{eV}^2$\t& "
	    << std::setw(4) << m12 / 1e-5 << " \\\\\n";
	out << "\t$\\sin^2 2\\theta_{12}$\t& "
	    << std::setw(4) << 4.0 * (1 - s12) * s12 << " \\\\\n";
	out << "\t\\midrule\n";
	out << "\t$\\Delta m_{23}^2/10^{-3}\\,\\text{eV}^2$\t& "
	    << std::setw(4) << m23 / 1e-3  << " \\\\\n";
	out << "\t$\\sin^2 2\\theta_{13}$\t& "
	    << std::setw(4) << 4.0 * (1 - s13) * s13 << " \\\\\n";
	out << "\t$\\sin^2 \\theta_{23}$\t& "
	    << std::setw(4) << s23 << " \\\\\n";
	out << "\t$\\delta_\\text{CP}$\t& "
	    << std::setw(4) << dCP << " \\\\\n";
	out << "\t\\bottomrule\n";
	out << "\\end{tabular}\n\n\n";

	out << "\\end{frame}\n\n";
}

void addTable(std::ofstream &out, std::string title, std::map<std::string, double> &mInt)
{
	out << "\\begin{frame}\n";
	addTable(out, title, "E", mInt);
	out << "\n\\vfill\n\n";
	addTable(out, title, "M", mInt);
	out << "\\end{frame}" << "\n\n";
}

void addTable(std::ofstream &out, std::string title, std::string p, std::map<std::string, double> &mInt)
{
	out << "\n{\\usebeamerfont{frametitle}\\usebeamercolor[fg]{frametitle} "
	    << "1R " << p << "-like, " << title << " mode}\n";
	out << "\\begin{center}\n";
	out << "\\footnotesize\n";

	out << std::setfill(' ') << std::setprecision(5) << std::fixed;
	out << "\\begin{tabular}{lrrrrrrr}" << "\n";
	out << "\t\\toprule" << "\n";
	out << "\t& $\\nu_\\mu$"
	     << "\t& $\\nu_e$"
	     << "\t& $\\overline{\\nu}_\\mu$"
	     << "\t& $\\overline{\\nu}_e$"
	     << "\t& $\\nu_\\mu \\to \\nu_e$"
	     << "\t& $\\overline{\\nu}_\\mu \\to \\overline{\\nu}_e$"
	     << "\t& Total"
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";
	//CCQE numbers
	out << "\tCCQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_"+p+"_CCQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCQE"] + mInt["nuE0_nuE0_"+p+"_CCQE"] +
				mInt["nuMB_nuMB_"+p+"_CCQE"] + mInt["nuEB_nuEB_"+p+"_CCQE"] +
				mInt["nuM0_nuE0_"+p+"_CCQE"] + mInt["nuMB_nuEB_"+p+"_CCQE"]
	     << "\t\\\\" << "\n";
	//CCnQE numbers
	out << "\tCCnQE\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_"+p+"_CCnQE"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCnQE"] + mInt["nuE0_nuE0_"+p+"_CCnQE"] +
				mInt["nuMB_nuMB_"+p+"_CCnQE"] + mInt["nuEB_nuEB_"+p+"_CCnQE"] +
				mInt["nuM0_nuE0_"+p+"_CCnQE"] + mInt["nuMB_nuEB_"+p+"_CCnQE"]
	     << "\t\\\\" << "\n";
	//NC numbers
	out << "\tNC\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_NC"] + mInt["nuE0_nuE0_"+p+"_NC"] +
				mInt["nuMB_nuMB_"+p+"_NC"] + mInt["nuEB_nuEB_"+p+"_NC"] +
				mInt["nuM0_nuE0_"+p+"_NC"] + mInt["nuMB_nuEB_"+p+"_NC"]
	     << "\t\\\\" << "\n";
	out << "\t\\midrule" << "\n";	
	//Total numbers
	out << "\tTot\t& "
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCQE"] + mInt["nuM0_nuM0_"+p+"_CCnQE"] + mInt["nuM0_nuM0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuE0_nuE0_"+p+"_CCQE"] + mInt["nuE0_nuE0_"+p+"_CCnQE"] + mInt["nuE0_nuE0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuMB_"+p+"_CCQE"] + mInt["nuMB_nuMB_"+p+"_CCnQE"] + mInt["nuMB_nuMB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuEB_nuEB_"+p+"_CCQE"] + mInt["nuEB_nuEB_"+p+"_CCnQE"] + mInt["nuEB_nuEB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuE0_"+p+"_CCQE"] + mInt["nuM0_nuE0_"+p+"_CCnQE"] + mInt["nuM0_nuE0_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuMB_nuEB_"+p+"_CCQE"] + mInt["nuMB_nuEB_"+p+"_CCnQE"] + mInt["nuMB_nuEB_"+p+"_NC"] << "\t&"
	     << std::setw(8) << mInt["nuM0_nuM0_"+p+"_CCQE"] + mInt["nuM0_nuM0_"+p+"_CCnQE"] + mInt["nuM0_nuM0_"+p+"_NC"] +
	     			mInt["nuE0_nuE0_"+p+"_CCQE"] + mInt["nuE0_nuE0_"+p+"_CCnQE"] + mInt["nuE0_nuE0_"+p+"_NC"] +
	     			mInt["nuMB_nuMB_"+p+"_CCQE"] + mInt["nuMB_nuMB_"+p+"_CCnQE"] + mInt["nuMB_nuMB_"+p+"_NC"] +
	     			mInt["nuEB_nuEB_"+p+"_CCQE"] + mInt["nuEB_nuEB_"+p+"_CCnQE"] + mInt["nuEB_nuEB_"+p+"_NC"] +
	     			mInt["nuM0_nuE0_"+p+"_CCQE"] + mInt["nuM0_nuE0_"+p+"_CCnQE"] + mInt["nuM0_nuE0_"+p+"_NC"] +
	     			mInt["nuMB_nuEB_"+p+"_CCQE"] + mInt["nuMB_nuEB_"+p+"_CCnQE"] + mInt["nuMB_nuEB_"+p+"_NC"] 
	     << "\t\\\\" << "\n";
	out << "\t\\bottomrule" << "\n";	
	out << "\\end{tabular}" << "\n";
	out << "\\end{center}" << "\n\n";
}
