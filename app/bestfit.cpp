#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

void beginDocument(std::ofstream &out);
void endDocument(std::ofstream &out);
void addFrame(std::ofstream &out, std::map<std::string, std::string> &plots);
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

	std::string fitInput;
	double fitPoint;

	cd->Get("fit_input", fitInput);
	cd->Get("fit_point", fitPoint);
	std::cout << "Files in " << fitInput << " for point " << fitPoint << std::endl;
	std::string cmd = "ls " + fitInput + " > .tmp_bestfit";
	std::cout << cmd << std::endl;
	system(cmd.c_str());

	int pp = fitPoint;

	double epsi[1000];
	int    Point;
	int    TPoint;

	bool keep = true;
	std::string file;
	std::ifstream listExclusion(".tmp_bestfit");
	while (keep && std::getline(listExclusion, file))
	{
		TFile *inf = new TFile(file.c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));

		sX2->SetBranchAddress("Epsilons", epsi);
		sX2->SetBranchAddress("Point",  &Point);
		sX2->SetBranchAddress("TPoint", &TPoint);

		for (int i = 0; i < sX2->GetEntries(); ++i)
		{
			sX2->GetEntry(i);

			if (Point == pp)
			{
				keep = false;
				break;
			}
		}
		inf->Close();
		delete inf;
	}

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
	//std::cout << osc->PMNS().cwiseAbs2() << std::endl;

	//std::string oscfile;
	//cd->Get("oscillation", oscfile);
	//TFile *oscr = new TFile(oscfile.c_str());

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

	//std::map<std::string, double> numevts, numsyst;
	Reco* reco;
	//TFile *rout;

	//TH1D* FCH1Re, *FHC1Rmu, *RHC1Re, *RHC1Rmu;
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

					std::cout << "adding to " << hname << std::endl;
					//std::string nn;
					//if (hname.find("_E_") != std::string::npos)
					//	nn = "E_" + horn[ih];
					//else if (hname.find("_M_") != std::string::npos)
					//	nn = "M_" + horn[ih];
					
					//if (hsys.count(nn))
					//	py->Multiply(hsys[nn]);

					if (refSpectra.count(hname) && refSpectra[hname])
						refSpectra[hname]->Add(py);
					else
						refSpectra[hname] = static_cast<TH1D*>(py->Clone());
					//numevts[hname] = py->Integral();

					//rout->cd();
					//py->SetName(hname.c_str());
					//py->Write(hname.c_str(), TObject::kWriteDelete);
				}
			}

			delete flux;
			delete reco;
		}

		//std::ofstream tout(name[ih].c_str());
		//tout << "\%" << horn[ih] << " mode\n\n";
		//output(tout, M12, M23, S12, S13, S23, dCP);
		//output(tout, numevts);
		////output(tout, numsyst);

		//tout.close();
		//rout->Close();

		//numevts.clear();
	}

	//systematics
	//std::map<std::string, std::string> syserre;
	//cd->Get("systematic_", syserre);
	//std::map<std::string, TH1D*> hsys;
	//std::cout << "applying syserre " << syserre.size() << std::endl;

	std::string baseOut, corrMatrix;
	cd->Get("base", baseOut);

	//load systematics in map
	std::map<std::string, TFile*> sysFile;
	std::map<std::string, TH1D*>::iterator it;
	for (it = refSpectra.begin(); it != refSpectra.end(); ++it)
	{
		std::string key = "systematic_" + it->first;
		std::string fileName;
		if (!cd->Get(key, fileName))
			continue;
		sysFile[it->first] = new TFile(fileName.c_str());
	}

	std::string texName = baseOut + "sys_variation.tex";
	std::ofstream tout(texName.c_str());
	beginDocument(tout);

	//get list of keys from first file, the name of systematics is the same
	TIter next(sysFile.begin()->second->GetListOfKeys());
	TKey *k;
	bool isSpline = false;
	int sN = 0;
	while ( k = static_cast<TKey*>(next()) )
	{
		std::map<std::string, std::string> plots;
		std::string sysname = k->GetName();
		std::cout << "systematic " << sysname << std::endl;

		if (sysname.find_first_of('_') != sysname.find_last_of('_'))	//spline systematic
		{
			std::cout << "spline" << std::endl;
			if (sysname.find("p1") != std::string::npos)
				isSpline = true;
			else
				continue;

			isSpline = true;
		}
		else
		{
			std::cout << "linear" << std::endl;
			isSpline = false;
		}


		if (!sysname.empty())
		{
			std::string hist_name = sysname;
			if (isSpline)
				sysname.erase(sysname.find_last_of('_'));

			for (it = refSpectra.begin(); it != refSpectra.end(); ++it)
			{
				//collect histogram with systematics
				//TH1D* hsys = static_cast<TH1D*>(sysFile[it->first]->Get(k->GetName()));
				//make a copy of histogram
				//sysSpectra[it->first] = static_cast<TH1D*>(it->second->Clone());

				TH1D* hsys = static_cast<TH1D*>(sysFile[it->first]->Get(hist_name.c_str()));
				//std::map<int, TH1D*> hsys;
				//for (is = spline.begin(); is != spline.end(); ++is)
				//	hsys[is->first] = static_cast<TH1D*>(sysFile[it->first]->Get(is->second.c_str()));


				int sysnum = std::stoi(sysname.substr(sysname.find_last_of('_')+1));
				double en = epsi[sysnum+67];
				std::string name = "bf_" + it->first + "_" + sysname;
				std::ofstream fout(name.c_str());
				fout << "# " << sysname << ", " << k->GetTitle() << " -> " << en << std::endl;
				std::cout << "look:  " << name << ", " << k->GetTitle() << " -> " << en << std::endl;

				//1, e, sigmas
				for (int b = 1; b < it->second->GetNbinsX()+1; ++b)
				{
					fout << it->second->GetBinLowEdge(b) << "\t" << it->second->GetBinContent(b);
					fout << "\t" << hsys->GetBinContent(b) * it->second->GetBinContent(b);
					fout << "\t" << 1 + en * (hsys->GetBinContent(b) - 1) * it->second->GetBinContent(b);
					fout << "\t" << hsys->GetBinContent(b);
					fout << "\t" << 1 + en * (hsys->GetBinContent(b) - 1);
					fout << std::endl;
				}

				fout.close();

				std::stringstream up;
				if (it->first.find("E") != std::string::npos)
					up << 1.2;
				else if (it->first.find("M") != std::string::npos)
					up << 10.0;

				std::string cmd = "gnuplot -e \'dataset=\"" + name +
						  "\"\' -e \'up=" + up.str() + "\' best.glp";
				std::cout << "Call " << cmd << std::endl;
				system(cmd.c_str());
				std::string move = "mv " + name + "* " + baseOut + "/";
				system(move.c_str());

				std::stringstream stt;
				stt << k->GetTitle();
				stt << ", $\\varepsilon = " << en << "$";
				plots[name] = stt.str();
				//std::string title = k->GetTitle();
			}

			addFrame(tout, plots);
		}
		++sN;
	}

	endDocument(tout);
	tout.close();

	chdir(baseOut.c_str());
	system("pdflatex sys_variation.tex");

	delete f1;

	return 0;
}

void beginDocument(std::ofstream &out)
{
	out << "\\documentclass[aspectratio=169]{beamer}\n";
	out << "\\setbeamertemplate{navigation symbols}{}\n";
	out << "\\setbeamerfont{frametitle}{size=\\small}\n";
	out <<"\n\\begin{document}\n\n";
}

void addFrame(std::ofstream &out, std::map<std::string, std::string> &plots)
{
	std::string title = plots.begin()->second;
	while (title.find('_') != std::string::npos)
		title.replace(title.find('_'), 1, " ");

	out << "\\begin{frame}\n";
	out << "\t\\frametitle{\\small " << title << "}\n";
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
		out << "\t\t\t\\resizebox{\\linewidth}{!}{\\input{" << name << ".tex}}\n";

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

void endDocument(std::ofstream &out)
{
	out <<"\n\\end{document}\n\n";
}
