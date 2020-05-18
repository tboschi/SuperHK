#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "event/Flux.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TMatrixD.h"
#include "TMatrixT.h"

#include "Eigen/Dense"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	std::string input;
	int truePoint;
	std::vector<int> fitPoints;

	cd->Get("input", input);
	cd->Get("true_point", truePoint);
	cd->Get("fit_point", fitPoints);
	fitPoints.insert(fitPoints.begin(), truePoint);	//truePoint at front

	int fitSize = fitPoints.size();

	std::cout << "Files in " << input << " for points " << truePoint << std::endl;
	std::string cmd = "ls " + input + " > .tmp_runchi2";
	system(cmd.c_str());

	int Point;
	double M12;
	double M23;
	double S12;
	double S13;
	double S23;
	double dCP;

	std::map<int, std::vector<double> >::iterator iparm;
	std::map<int, std::vector<double> > tparm;

	std::string file;
	std::ifstream listExclusion(".tmp_runchi2");
	//looking for fit points
	std::vector<int> findFits = fitPoints; //make a copy to find fit points
	while (findFits.size() && std::getline(listExclusion, file))
	{
		TFile *inf = new TFile(file.c_str(), "READ");
		if (inf->IsZombie())
		{
			std::cout << "file doesn't exist\n";
			continue;
		}

		TTree* sX2 = static_cast<TTree*>(inf->Get("mcTree"));

		sX2->SetBranchAddress("Point", &Point);
		sX2->SetBranchAddress("M12", &M12);
		sX2->SetBranchAddress("M23", &M23);
		sX2->SetBranchAddress("S12", &S12);
		sX2->SetBranchAddress("S13", &S13);
		sX2->SetBranchAddress("S23", &S23);
		sX2->SetBranchAddress("CP",  &dCP);

		sX2->GetEntry(0);
		int p0 = Point;				//first point
		sX2->GetEntry(sX2->GetEntries()-1);
		int p1 = Point;				//last point

		//check any point is in this file
		int jp = -1;
		std::vector<int>::iterator ip;
		for (ip = findFits.begin(); ip != findFits.end(); )
		{
			if (p0 <= *ip && *ip <= p1)
			{
				jp = *ip;
				ip = findFits.erase(ip);
				break;
			}
			else
				++ip;
		}

		if (jp >= 0)	//if there is a point loop and find parameters
			for (int i = 0; i < sX2->GetEntries(); ++i)
			{
				sX2->GetEntry(i);

				if (jp == Point && !tparm.count(Point))
				{
					std::cout << "Found point " << Point << std::endl;

					std::vector<double> tv;
					tv.push_back(M12);	//entry 0
					tv.push_back(M23);	//entry 1
					tv.push_back(S12);	//entry 2
					tv.push_back(S13);	//entry 3
					tv.push_back(S23);	//entry 4
					tv.push_back(dCP);	//entry 5

					tparm[Point] = tv;
				}
			}

		inf->Close();
		delete inf;
	}


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

	std::map<std::string, TH1D*> allSpectra;

	//looping on each fit point, creating the erco spectra and saving them to allSpectra
	for (iparm = tparm.begin(); iparm != tparm.end(); ++iparm)
	{
		std::vector<double> tv = iparm->second;

		Oscillator *osc = new Oscillator(densityFile);
		osc->SetMasses<Oscillator::normal>(tv.at(0), tv.at(1));
		osc->SetPMNS<Oscillator::sin2>(tv.at(2), tv.at(3), tv.at(4), tv.at(5));

		for (int ih = 0; ih < 2; ++ih)	//FHC, RHC
			for (int im = 0; im < 6; ++im)	//nuE->nuE, nuM->nuE, nuM->nuM
		{
			cd->Get("reco_" + mode[im] + "_" + horn[ih], reco_file);
			Reco *reco = new Reco(reco_file);

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


				TH1D *py = reco->Project(chan[ic], 'y');
				if (py)
				{
					std::string hname = chan[ic].substr(0, chan[ic].find_first_of('_'))
							    + "_" + horn[ih] + "_" + std::to_string(iparm->first);

					//std::cout << "adding to " << hname << std::endl;
					//std::string nn;
					//if (hname.find("_E_") != std::string::npos)
					//	nn = "E_" + horn[ih];
					//else if (hname.find("_M_") != std::string::npos)
					//	nn = "M_" + horn[ih];

					//if (hsys.count(nn))
					//	py->Multiply(hsys[nn]);

					if (allSpectra.count(hname) && allSpectra[hname])
						allSpectra[hname]->Add(py);
					else
						allSpectra[hname] = static_cast<TH1D*>(py->Clone());
					//numevts[hname] = py->Integral();

					//py->SetName(hname.c_str());
					//py->Write(hname.c_str(), TObject::kWriteDelete);
				}
			}

			delete flux;
			delete reco;
		}
	}

	//systematics
	//std::map<std::string, std::string> syserre;
	//cd->Get("systematic_", syserre);
	//std::map<std::string, TH1D*> hsys;
	//std::cout << "applying syserre " << syserre.size() << std::endl;

	//load systematics in map with key E_FHC, E_RHC, etc..
	std::map<std::string, std::vector<TH1D*> > sysMap;
	std::map<std::string, TH1D*>::iterator it;
	for (it = allSpectra.begin(); it != allSpectra.end(); ++it)
	{
		std::string stype = it->first.substr(0, it->first.find_last_of('_'));
		if (sysMap.count(stype))
				continue;

		std::string fileName;
		if (!cd->Get("systematic_" + stype, fileName))
			continue;

		TFile * sysFile = new TFile(fileName.c_str());
		TIter next(sysFile->GetListOfKeys());
		TKey *k;
		std::string spline;
		//bool isSpline = false;
		std::vector<TH1D*> sysHistograms;
		while ( k = static_cast<TKey*>(next()) )
		{
			std::string sysname = k->GetName();
			std::cout << "systematic " << sysname << std::endl;

			if (sysname.find_first_of('_') != sysname.find_last_of('_'))	//spline systematic
			{
				if (sysname.find("p1") != std::string::npos)
					spline = sysname;
			}
			else
				spline = sysname;

			if (spline.empty())
				continue;

			TH1D* hsys = static_cast<TH1D*>(sysFile->Get(spline.c_str()));
			hsys->SetDirectory(0);
			sysHistograms.push_back(hsys);
		}

		sysMap[stype] = sysHistograms;
		sysFile->Close();
	}

	std::cout << "found systematics" << std::endl;

	std::string matrixFile;
	cd->Get("correlation", matrixFile);
	TFile * mF = new TFile(matrixFile.c_str());
	TMatrixT<double> * cmat = static_cast<TMatrixT<double>*>(mF->Get("correlation"));

	Eigen::MatrixXd corr(cmat->GetNcols(), cmat->GetNcols());
	for (int r = 0; r < cmat->GetNcols(); ++r)
		for (int c = 0; c < cmat->GetNcols(); ++c)
		corr(r, c) = cmat->operator()(r, c);

	std::string baseName = std::to_string(fitPoints.front()); //has a trailing _
	for (int t = 1; t < fitPoints.size(); ++t)
		baseName += "_" + std::to_string(fitPoints[t]);

	std::string base;
	cd->Get("base", base);
	std::ofstream tout(base.c_str());

	for (double e1 = -10; e1 < 10.5; e1 += 0.5)
		for (double e2 = -10; e2 < 10.5; e2 += 0.5)
	{
		double yn = 0;

		Eigen::VectorXd sysValues(corr.cols());
		sysValues[0] = e1;
		sysValues[1] = e2;

		//double chi2_0 = sysValues.transpose() * corr.inverse() * sysValues;
		double chi2_0 = sysValues.transpose() * sysValues;
		std::vector<double> chi2(fitPoints.size() - 1, chi2_0);

		std::map<std::string, std::vector<TH1D*> >::iterator is;
		for (is = sysMap.begin(); is != sysMap.end(); ++is)
		{
			TH1D* hObs = allSpectra[is->first + "_" + std::to_string(fitPoints.front())];
			for (int b = 1; b < hObs->GetNbinsX()+1; ++b)
			{
				double On = hObs->GetBinContent(b);

				double yn = 0;
				for (int i = 0; i < sysValues.size(); ++i)
				{
					TH1D* hsys = (is->second)[i];
					yn += sysValues[i] * (hsys->GetBinContent(b) - 1);
				}

				for (int t = 1; t < fitPoints.size(); ++t)
				{
					TH1D* hExp = allSpectra[is->first + "_" + std::to_string(fitPoints[t])];
					double En = hExp->GetBinContent(b);

					double chi2n = 2 * (En * (1 + yn) - On);
					if (On > 0 && En > 0)
						chi2n += 2 * On * log(On / (1 + yn) / En);

					chi2[t-1] += chi2n;
				}
			}
		}

		tout << e1 << "\t" << e2;
		for (int t = 1; t < fitPoints.size(); ++t)
			tout << "\t" << chi2[t-1];
		tout << std::endl;
	}


	tout.close();
	delete f1;

	return 0;
}
