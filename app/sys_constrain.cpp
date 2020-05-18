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

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	//get systematics and correlation matrix

	std::string baseOut, corrMatrix;
	cd->Get("base", baseOut);

	int nBin = NumBin(cd);
	int nSys = NumSys(cd);

	//load systematics in matrix
	std::string type[4] = {"E_FHC", "E_RHC", "M_FHC", "M_RHC"};

	Eigen::MatrixXd sysMatrix = Eigen::MatrixXd::Zero(4 * nBin, nSys);
	std::map<int, Eigen::MatrixXd> splineMatrices;
	for (t = 0; t < 4; ++t)
	{
		std::string fileName;
		if (!cd->Get("systematic_" + type[t], fileName))
			continue;

		TFile *sysf = new TFile(fileName.c_str());
		TIter next(sysf->GetListOfKeys());
		TKey *k;
		std::map<int, std::string> spline;
		std::map<int, std::string>::iterator ip;
		while (k = static_cast<TKey*>(next()))
		{
			std::string sysname = k->GetName();
			std::cout << "accepting " << sysname << std::endl;

			if (sysname.find_first_of('_') != sysname.find_last_of('_'))	//spline systematic
			{
				if (sysname.find("m3") != std::string::npos)
					sigma = -3;
				else if (sysname.find("m1") != std::string::npos)
					sigma = -1;
				else if (sysname.find("p1") != std::string::npos)
					sigma = 1;
				else if (sysname.find("p3") != std::string::npos)
					sigma = 3;

				spline[sigma] = sysname;
			}
			else
				spline[1] = sysname;

			if ( (spline.size() == 4) || (spline.size() == 1) )
			{		//else do nothing

				Eigen::MatrixXd syst = Eigen::MatrixXd::Zero(4 * nBin, spline.size());
				for (ip = spline.begin(); ip != spline.end(); ++ip)
				{
  					TH1D* hsys = static_cast<TH1D*>(sysf->Get(ip->second.c_str()));
					for (int n = 0; n < hsys->GetNbinsX(); ++n)
						syst(nBin * t + n, ip->first) = (hsys->GetBinContent(n+1) - 1);
				}

				sysname = spline.begin()->second;
				if (sysname.find_first_of('_') != sysname.find_last_of('_'))	//spline
					sysname.erase(sysname.find_last_of('_'));
				sysname.erase(0, sysname.find_first_of('_'));	//sysname is just number

				//number of systematic	
				int k = std::strtol(sysname, NULL, 10);

				int p1col = spline.size() == 1 ? 0 : 2;

				sysMatrix.col(k) += syst.col(p1col);	//add 1sigma to main matrix

				if (spline.size() > 1)
				{	//then add spline matrices to map
					if (splineMatrices.count(k))
						splineMatrice[k] += syst;
					else
						splineMatrice[k] = Eigen::MatrixXd::Zero(4 * nBin, spline.size());
				}

				spline.clear();
			}
		}
	}


	//systematics loaded,
	//now move on to creating spectra and looping over files

	double M12;
	double M23;
	double S12;
	double S13;
	double S23;
	double dCP;

	double epsil[200];
	double X2;	//min value of X2 at given point

	std::string file;
	std::ifstream listExclusion(".tmp_constrain");

	//get first file for True point
	std::getline(listExclusion, file);
	TFile *inf = new TFile(file.c_str(), "READ");
	if (inf->IsZombie())
	{
		std::cout << "file doesn't exist\n";
		continue;
	}

	TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));

	sX2->SetBranchAddress("TM12", &M12);
	sX2->SetBranchAddress("TM23", &M23);
	sX2->SetBranchAddress("TS12", &S12);
	sX2->SetBranchAddress("TS13", &S13);
	sX2->SetBranchAddress("TS23", &S23);
	sX2->SetBranchAddress("TCP",  &dCP);

	sX2->GetEntry(0);

	Oscillator *osc = new Oscillator(densityFile);

	osc->SetMasses<Oscillator::normal>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	//spectra for true point ( O_n )
	//std::map<std::string, TH1D*> trueSpectra = ConstructSpectrum(cd, osc);
	//it is a concat of all events, E_FHC + E_RHC + M_FHC + M_RHC
	Eigen::VectorXd trueSpectra = ConstructSpectrum(cd, osc);

	inf->Close();

	//now loop over all files to get the fit spectra ( E_n )
	while (std::getline(listExclusion, file))
	{
		inf = new TFile(file.c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* sX2 = static_cast<TTree*>(inf->Get("stepX2Tree"));

		sX2->SetBranchAddress("Espilons", epsil);
		sX2->SetBranchAddress("X2", &X2);

		sX2->SetBranchAddress("M12", &M12);
		sX2->SetBranchAddress("M23", &M23);
		sX2->SetBranchAddress("S12", &S12);
		sX2->SetBranchAddress("S13", &S13);
		sX2->SetBranchAddress("S23", &S23);
		sX2->SetBranchAddress("CP",  &dCP);

		osc->SetMasses<Oscillator::normal>(M12, M23);
		osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

		for (int i = 0; i < sX2->GetEntries(); ++i)
		{
			sX2->GetEntry(i);

			//std::map<std::string, TH1D*> fitSpectra = ConstructSpectrum(cd, osc);
			Eigen::VectorXd fitSpectra = ConstructSpectrum(cd, osc);

			// compute chi2
			for (int n = 0; n < 4 * nBin; ++n)
			{
				double On = trueSpectra(n);

				double sigma = 1;
				//remove unit so it can be multiplied by sigmas
				double yn = sigma * (hsys->GetBinContent(b) - 1);

				fout << hObs->GetBinLowEdge(b) << "\t" << On;
				for (int t = 1; t < fitPoints.size(); ++t)
				{
					//these are the expected events (= fit points)
					TH1D* hExp = allSpectra[is->first + "_" + std::to_string(fitPoints[t])];
					double En = hExp->GetBinContent(b);
					double chi = 0;

					double chi2 = 2 * (En * (1 + yn) - On);
					if (On > 0 && En > 0)
						chi2 += 2 * On * log(On / (1 + yn) / En);
					totChi2[t-1] += chi2;

					fout << "\t" << En << "\t" << chi2;
				}
				fout << std::endl;
			}

		}

		inf->Close();
	}


		//loop over fij hist / or loop over each file and get the same syserre
		std::map<std::string, TFile*>::iterator is;
		for (is = sysFile.begin(); is != sysFile.end(); ++is)
		{
			TH1D* hsys = static_cast<TH1D*>(is->second->Get(spline.c_str()));

			std::string name = baseName + "_" + is->first + "_" + sysname;
			std::ofstream fout(name.c_str());
			std::cout << "look:  " << name << ", " << k->GetTitle() << std::endl;
			fout << "# " << sysname << ", " << k->GetTitle() << std::endl;
			fout << "#bin\t" << std::to_string(fitPoints.front());
			for (int t = 1; t < fitPoints.size(); ++t)
				fout << "\t" << fitPoints[t];
			fout << std::endl;

			//these are the observed events (= true point)
			TH1D* hObs = allSpectra[is->first + "_" + std::to_string(fitPoints.front())];
			std::vector<double> totChi2(fitPoints.size() - 1);
			for (int b = 1; b < hObs->GetNbinsX()+1; ++b)
			{
				double On = hObs->GetBinContent(b);

				double sigma = 1;
				//remove unit so it can be multiplied by sigmas
				double yn = sigma * (hsys->GetBinContent(b) - 1);

				fout << hObs->GetBinLowEdge(b) << "\t" << On;
				for (int t = 1; t < fitPoints.size(); ++t)
				{
					//these are the expected events (= fit points)
					TH1D* hExp = allSpectra[is->first + "_" + std::to_string(fitPoints[t])];
					double En = hExp->GetBinContent(b);
					double chi = 0;

					double chi2 = 2 * (En * (1 + yn) - On);
					if (On > 0 && En > 0)
						chi2 += 2 * On * log(On / (1 + yn) / En);
					totChi2[t-1] += chi2;

					fout << "\t" << En << "\t" << chi2;
				}
				fout << std::endl;
			}

			fout << "\n#tot chi2";
			for (int c = 0; c < totChi2.size(); ++c)
				fout << "\t" << totChi2[c];
			fout << std::endl;

			fout.close();

			std::stringstream up;
			if (is->first.find("E") != std::string::npos)
				up << 1.25;
			else if (is->first.find("M") != std::string::npos)
				up << 10.0;

			cmd = "gnuplot -e \'dataset=\"" + name +
				"\"\' -e \'up=" + up.str() + "\' chi2.glp";
			std::cout << "Call " << cmd << std::endl;
			//system(cmd.c_str());
			std::string move = "mv " + name + "* " + baseOut + "/";
			system(move.c_str());

			plots[name] = k->GetTitle();
			//std::string title = k->GetTitle();
		}

		spline.clear();
		addFrame(tout, plots);
	}

	endDocument(tout);
	tout.close();

	cmd = "pdflatex " + baseName + "_chi2_variation.tex";
	chdir(baseOut.c_str());
	//system(cmd.c_str());

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

int NumBins(CardDealer *cd)
{
	std::string reco_file;

	if (!cd->Get("reco_nuE0_nuE0_FHC", reco_file))
		return 0;

	Reco *reco = new Reco(reco_file);
	const double *bins;
	int nBin = reco->BinsX(bins);

	delete reco;

	return nBin;
}

int NumSys(CardDealer *cd)
{
	std::string fileName;
	if (!cd->Get("systematic_" + stype, fileName))
		return 0;

	TFile *sysf = new TFile(fileName.c_str());
	TIter next(sysf->GetListOfKeys());
	TKey *k;
	std::string sysname_;
	int count = 0;
	while (k = static_cast<TKey*>(next()))
	{
		sysname = k->GetName();
		if (sysname.find_first_of('_') != sysname.find_last_of('_'))	//spline systematic
			sysname.erase(sysname.find_last_of('_'));

		if (sysname.find(sysname_) == std::string::npos) //systematic different from previous
			++count;

		sysname_ = sysname;
	}

	sysf->Close();

	return count;
}

//std::map<std::string, TH1D*> ConstructSpectrum(CardDealer *cd, Oscillator *osc)
Eigen::VectorXd ConstructSpectrum(CardDealer *cd, Oscillator *osc)
{
	std::string reco_file;

	std::string mode[6] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string chan[6] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	std::string horn[2] = {"FHC", "RHC"};

	Nu fin[6]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu fout[6] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};
	std::string oscf[6] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	std::map<std::string, TH1D*> spectra;
	std::map<std::string, TH1D*>::iterator is;

	int nBin = 0;
	for (int ih = 0; ih < 2; ++ih)	//FHC, RHC
		for (int im = 0; im < 6; ++im)	//nuE->nuE, nuM->nuE, nuM->nuM
	{
		cd->Get("reco_" + mode[im] + "_" + horn[ih], reco_file);
		Reco *reco = new Reco(reco_file);

		const double *bins;
		nBin = reco->BinsX(bins);
		f1->SetRange(bins[0], bins[nBin-1]);

		std::string fname = horn[ih] + "_" + mode[im];
		TH1D* flux = new TH1D(fname.c_str(), fname.c_str(), nBin, bins);
		//TH1D* flux = static_cast<TH1D*>(oscr->Get(oscf[im].c_str()));
		flux->SetDirectory(0);
		flux->Add(f1);	//filled with ones

		osc->Oscillate(fin[im], fout[im], flux);

		for (int ic = 0; ic < 6; ++ic)	//CCQE, CCnQE, NC x E, M
		{
			if (chan[ic].find("NC") == std::string::npos)	//it is not a NC
				reco->Scale(chan[ic], flux);
			else if (mode[im] == "nuM0_nuE0" || mode[im] == "nuMB_nuEB") //and no events for
				reco->Scale(chan[ic], 0.0);			     //these two
			//else no need to scale for NC event
			//reco->Scale(chan[ic], flux);

			TH1D *py = reco->Project(chan[ic], 'y');
			if (py)
			{
				std::string hname = chan[ic].substr(0, chan[ic].find_first_of('_'))
					+ "_" + horn[ih] + "_" + std::to_string(iparm->first);

				if (spectra.count(hname) && spectra[hname])
					spectra[hname]->Add(py);
				else
					spectra[hname] = static_cast<TH1D*>(py->Clone());
			}
		}

		delete flux;
		delete reco;
	}

	delete f1;

	std::string type[4] = {"E_FHC", "E_RHC", "M_FHC", "M_RHC"};

	Eigen::VectorXd vect(4 * nBin);
	nBin = 0;
	for (t = 0; t < 4; ++t)
	{
		for (int i = 0; i < spectra[type[t]]->GetNbinsX(); ++i) 
			vect(i-1) = spectra[type[t]]->GetBinContent(i+1);

		totBins += spectra[type[t]]->GetNbinsX();

		delete spectra[type[t]];
	}

	return vect;
}
