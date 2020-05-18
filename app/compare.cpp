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

void addTable(std::ofstream &out, std::string type,
	      std::map<std::string, double> &mInt,
	      std::vector<double> &parms);
void addTable(std::ofstream &out, std::string type, std::string p,
	      std::map<std::string, double> &mInt);
void addFrame(std::ofstream &out, std::vector<double> &parms,
	      std::map<std::string, std::string> &plots);
int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	TRandom3 *mt = new TRandom3();
	CardDealer *cd = new CardDealer(cardFile);

	//std::string texFHC, texRHC;
	//cd->Get("out_0", texFHC);
	//cd->Get("out_B", texRHC);
	//texFHC.erase(texFHC.find(".root"));
	//texRHC.erase(texRHC.find(".root"));

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	double M12, M23;
	double S12, S13, S23;
	double dCP;
	std::string s_M12, s_M23;
	std::string s_S12, s_S13, s_S23;
	std::string s_dCP;

	cd->Get("M12", s_M12);
	cd->Get("M23", s_M23);
	cd->Get("S12", s_S12);
	cd->Get("S13", s_S13);
	cd->Get("S23", s_S23);
	cd->Get("dCP", s_dCP);

	s_M12.erase(0, 1);
	s_M23.erase(0, 1);
	s_S12.erase(0, 1);
	s_S13.erase(0, 1);
	s_S23.erase(0, 1);
	s_dCP.erase(0, 1);

	M12 = std::stod(s_M12);
	M23 = std::stod(s_M23);
	S12 = std::stod(s_S12);
	S13 = std::stod(s_S13);
	S23 = std::stod(s_S23);
	dCP = std::stod(s_dCP);

	std::vector<double> parms;
	parms.push_back(M12);
	parms.push_back(M23);
	parms.push_back(S12);
	parms.push_back(S13);
	parms.push_back(S23);
	parms.push_back(dCP);

	std::map<std::string, TH1D*> valorSpectra;
	std::string input;
	cd->Get("input", input);

	std::string file = input + "_dm21sq_"  + s_M12 +
				   "_dm32sq_"  + s_M23 +
				   "_theta12_" + s_S12 +
				   "_theta23_" + s_S23 +
				   "_cpv_"     + s_dCP +
				   "_theta13_" + s_S13 +
				   ".root";

	TFile *inf = new TFile(file.c_str(), "READ");

	if (!inf || inf->IsZombie())
	{
		std::cout << "file doesn't exist\n";
		return 1;
	}

	valorSpectra["E_FHC"] = static_cast<TH1D*> (inf->Get("valor_pdfts_ob_rf_superk__JPARC_FHC__1re_Ereco"));
	valorSpectra["E_RHC"] = static_cast<TH1D*> (inf->Get("valor_pdfts_ob_rf_superk__JPARC_RHC__1re_Ereco"));
	valorSpectra["M_FHC"] = static_cast<TH1D*> (inf->Get("valor_pdfts_ob_rf_superk__JPARC_FHC__1rmu_Ereco"));
	valorSpectra["M_RHC"] = static_cast<TH1D*> (inf->Get("valor_pdfts_ob_rf_superk__JPARC_RHC__1rmu_Ereco"));

	Oscillator *osc = new Oscillator(densityFile);
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

	if (file.find(".root") != std::string::npos)
		file.erase(file.find(".root"));
	if (file.find_last_of('/') != std::string::npos)
		file.erase(0, file.find_last_of('/')+1);
	while (file.find('.') != std::string::npos)
		file.replace(file.find('.'), 1, "_");

	std::string baseOut;
	cd->Get("base", baseOut);

	std::string texName = baseOut + file + ".tex";
	std::ofstream tout(texName.c_str());

	Reco* reco;

	std::map<std::string, double> numevts;
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
							    + "_" + horn[ih];

					std::string nname = mode[im] + "_" + chan[ic];
					numevts[nname] = py->Integral();

					if (refSpectra.count(hname) && refSpectra[hname])
						refSpectra[hname]->Add(py);
					else
						refSpectra[hname] = static_cast<TH1D*>(py->Clone());

				}
			}

			delete flux;
			delete reco;
		}

		addTable(tout, horn[ih], numevts, parms);
		numevts.clear();
	}


	std::map<std::string, std::string> plots;
	std::map<std::string, TH1D*>::iterator it;
	for (it = refSpectra.begin(); it != refSpectra.end(); ++it)
	{
		//collect histogram with systematics

		std::string name = file + "_" + it->first;
		std::string datName = name + ".dat";
		std::ofstream fout(datName.c_str());

		TH1D* hval = valorSpectra[it->first];
		for (int b = 1; b < it->second->GetNbinsX()+1; ++b)
			fout << it->second->GetBinLowEdge(b) << "\t"
			     << it->second->GetBinContent(b) << "\t"
			     << hval->GetBinContent(b) << std::endl;

		fout.close();

		std::stringstream up;
		if (it->first.find("E") != std::string::npos)
			up << 1.25;
		else if (it->first.find("M") != std::string::npos)
			up << 7.0;

		std::string cmd = "gnuplot -e \'dataset=\"" + name +
			"\"\' -e \'up=" + up.str() + "\' spectra.glp";
		//std::cout << "Call " << cmd << std::endl;
		system(cmd.c_str());
		std::string move = "mv " + name + "* " + baseOut + "/";
		system(move.c_str());

		plots[it->first] = name;
		//std::string title = k->GetTitle();
	}

	addFrame(tout, parms, plots);

	tout.close();
	inf->Close();

	std::cout << texName << std::endl;

	delete f1;

	return 0;
}


void addTable(std::ofstream &out, std::string type, std::map<std::string, double> &mInt, std::vector<double> &parms)
{
	std::stringstream title;
	title << "$\\Delta m_{12}^2/ 10^{-5}\\,\\text{eV}^2 = "
	      << std::setw(4) << parms[0] / 1e-5 << "$, ";
	title << "$\\sin^2 2\\theta_{12} = "
	      << std::setw(4) << 4.0 * (1 - parms[2]) * parms[2] << "$, ";
	title << "$\\Delta m_{23}^2/10^{-3}\\,\\text{eV}^2 = "
	      << std::setw(4) << parms[1] / 1e-3  << "$, ";
	title << "$\\sin^2 2\\theta_{13} = "
	      << std::setw(4) << 4.0 * (1 - parms[3]) * parms[3] << "$, ";
	title << "$\\sin^2 \\theta_{23} = "
	      << std::setw(4) << parms[4] << "$, ";
	title << "$\\delta_\\text{CP} = "
	      << std::setw(4) << parms[5] << "$";

	out << "\\begin{frame}\n";
	out << "\t\\frametitle{\\small " << title.str() << "}\n";
	out << "\t\\footnotesize\n";
	out << "\t\\centering\n";

	addTable(out, type, "E", mInt);

	out << "\n\\vfill\n\n";

	addTable(out, type, "M", mInt);

	out << "\\end{frame}" << "\n\n";
}

void addTable(std::ofstream &out, std::string type, std::string p, std::map<std::string, double> &mInt)
{
	out << "\n\\textbf{1R " << p << "-like, " << type << " mode}\n\n\\medskip\n";
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

void addFrame(std::ofstream &out, std::vector<double> &parms, std::map<std::string, std::string> &plots)
{
	out << "\% You must create a document that will input this file with\n";
	out << "\n\%\t\\input{filename}\n\n";
	out << "\% In the preamble put the following lines\n\n";
	out << "\%\t\\documentclass[aspectratio=169]{beamer}\n";
	out << "\%\t\\setbeamertemplate{navigation symbols}{}\n";
	out << "\%\t\\setbeamerfont{frametitle}{size=\\small}\n\n";
	out <<"\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\\\n\n";

	std::stringstream title;
	title << "$\\Delta m_{12}^2/ 10^{-5}\\,\\text{eV}^2 = "
	      << std::setw(4) << parms[0] / 1e-5 << "$, ";
	title << "$\\sin^2 2\\theta_{12} = "
	      << std::setw(4) << 4.0 * (1 - parms[2]) * parms[2] << "$, ";
	title << "$\\Delta m_{23}^2/10^{-3}\\,\\text{eV}^2 = "
	      << std::setw(4) << parms[1] / 1e-3  << "$, ";
	title << "$\\sin^2 2\\theta_{13} = "
	      << std::setw(4) << 4.0 * (1 - parms[3]) * parms[3] << "$, ";
	title << "$\\sin^2 \\theta_{23} = "
	      << std::setw(4) << parms[4] << "$, ";
	title << "$\\delta_\\text{CP} = "
	      << std::setw(4) << parms[5] << "$";

	out << "\\begin{frame}\n";
	out << "\t\\frametitle{\\small " << title.str() << "}\n";
	out << "\t\\footnotesize\n";
	out << "\t\\centering\n";

	out << "\t\\begin{columns}[c]\n";

	std::map<std::string, std::string>::iterator ip;
	int i = 0;
	for (ip = plots.begin(); ip != plots.end(); ++ip, ++i)
	{
		if (i % 2 == 0)
			out << "\t\t\\begin{column}{0.5\\linewidth}\n";

		std::string name = ip->first;
		if (name.find('/') != std::string::npos)
			name.erase(0, name.find_last_of('/')+1);
		std::string type = name;
		while (type.find('_') != std::string::npos)
			type.replace(type.find('_'), 1, " ");

		out << "\t\t\t\\textbf{" << type << "}\n";
		out << "\t\t\t\\resizebox{\\linewidth}{!}{\\input{" << ip->second << ".tex}}\n";

		if (i % 2 == 0)
			out << "\t\t\t\\vfill\n";
		else if (i % 2 == 1)
			out << "\t\t\\end{column}\n";

		if (i == 1)
			out << "\t\t\\hfill\n";
	}

	out << "\t\\end{columns}\n";
	out << "\\end{frame}\n\n";
}

