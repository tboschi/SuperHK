#include <fstream>
#include <iostream>
#include <iomanip>

#include <omp.h>

#include "tools/CardDealer.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

#include "TF1.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

void output(std::ofstream &out, std::map<std::string, double> &mInt);
int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	std::string treeName, files, outname, nameList;
	cd->Get("tree", treeName);
	cd->Get("input", files);
	cd->Get("output", outname);

	std::string cmd = "ls " + files + " > .tmp_list";
	system(cmd.c_str());
	std::ifstream in(".tmp_list");
	std::vector<std::string> fileList;
	while (std::getline(in, nameList))
		fileList.push_back(nameList);

	std::string mode[6] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string chan[6] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	std::string horn[2] = {"FHC", "RHC"};
	//std::string name[2] = {texFHC+".tex",  texRHC+".tex"};
	//std::string root[2] = {texFHC+".root", texRHC+".root"};
	//std::string oscf[6] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};
	int offs[6] = {2, 0, 2, 0, 2, 0};

	Nu fin[6]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu fout[6] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};

	std::string reco_file;
	std::vector<Reco*> vReco;
	for (int ih = 0; ih < 2; ++ih)
	{
		for (int im = 0; im < 6; ++im)
		{
			cd->Get("reco_" + mode[im] + "_" + horn[ih], reco_file);
			vReco.push_back(new Reco(reco_file));
		}
	}

	const double *bins;
	int nBinBeam = vReco.front()->BinsY(bins);
	int nBinEnu  = vReco.front()->BinsX(bins);
	//std::cout << "Creating " << nBinBeam << " bins" << std::endl;

	//TChain *ch = new TChain(treeName.c_str());
	//ch->Add(files.c_str());

	int nBinAtm = 2224;
	int nBinTot = nBinAtm + nBinBeam * 8;
	//  8 = 2 * 4
	double *mcBins = new double[nBinAtm + nBinBeam * 8];
	double *newmc  = new double[nBinAtm + nBinBeam * 8];

	int Point;
	double M12, M23;
	double S12, S13, S23;
	double dCP;

	std::string densityHK, densityKD;
	cd->Get("densityHK", densityHK);
	cd->Get("densityKD", densityKD);

	std::string MH;
	cd->Get("hierarchy", MH);

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);
	f1->SetRange(bins[0], bins[nBinBeam-1]);

	//std::cout << "allocating " << ssl.str() << std::endl;

	for (int f = 0; f < fileList.size(); ++f)
	{
		std::string name(fileList[f]);
		std::cout << "Reading " << name << "\t";
		std::string base = name.substr(name.find_last_of('/')+1);
		base = outname + "/" + base;
		std::cout << "into " << base << "\n";

		TFile input(name.c_str());

		TTree *oldt = static_cast<TTree*> (input.Get(treeName.c_str()));

		oldt->SetBranchAddress("mcBins", mcBins);
		oldt->SetBranchAddress("Point", &Point);
		oldt->SetBranchAddress("CP",  &dCP); 
		oldt->SetBranchAddress("M12", &M12);
		oldt->SetBranchAddress("M23", &M23);
		oldt->SetBranchAddress("S12", &S12);
		oldt->SetBranchAddress("S13", &S13);
		oldt->SetBranchAddress("S23", &S23);

		TFile output(base.c_str(), "RECREATE");
		output.cd();

		TTree *newt = new TTree("mcTree", oldt->GetTitle());

		std::stringstream ssl;
		ssl << nBinAtm + nBinBeam * 8;
		std::string mcbin = "mcBins[" + ssl.str() + "]/D";
		newt->Branch("mcBins", newmc, mcbin.c_str());
		newt->Branch("Point", &Point, "Point/I");
		newt->Branch("CP",  &dCP, "CP/D"); 
		newt->Branch("M12", &M12, "M12/D");
		newt->Branch("M23", &M23, "M23/D");
		newt->Branch("S12", &S12, "S12/D");
		newt->Branch("S13", &S13, "S13/D");
		newt->Branch("S23", &S23, "S23/D");

		double _dCP = -9999;
		std::cout << "Number of entries: " << oldt->GetEntriesFast() << std::endl;
		for (int n = 0; n < oldt->GetEntriesFast(); ++n)
		{
			std::cout << "Entry " << n << std::endl;
			oldt->GetEntry(n);

			if (dCP != _dCP)
			{
				//std::cout << "dCP --- " << dCP << std::endl;
				_dCP = dCP;
			}

			//copy atmospheric
			memcpy(newmc, mcBins,      sizeof(double) * nBinAtm);
			memset(newmc + nBinAtm, 0, sizeof(double) * nBinBeam * 8);

			//double number[4];
			//for (int i = 0; i < 4; ++i)
			//	number[i] = 0;
			//generate beam
			//
			std::vector<Reco*> copyReco;
			std::vector<TH1D*> vFlux;
			for (int jj = 0; jj < 2*vReco.size(); ++jj)
			{
				int ih = (jj / 2) / 6;
				int im = (jj / 2) % 6;
				std::string detName;
				if (jj % 2 == 0)
					detName = "HK";
				else
					detName = "KD";

				std::string fname = horn[ih] + "_" + mode[im] + "_" + detName;

				copyReco.push_back(new Reco(*vReco.at(jj/2)));
				vFlux.push_back(new TH1D(fname.c_str(), fname.c_str(), nBinEnu, bins));
			}

//#pragma omp parallel for num_threads(2) default(none) shared(std::cout, densityHK, densityKD, nBinAtm, nBinBeam, f1, copyReco, vFlux, dCP, M12, M23, S12, S13, S23, MH, mode, chan, offs, newmc, fout, fin)
			for (int jj = 0; jj < copyReco.size(); ++jj)
			{
				int ih = (jj / 2) / 6;
				int im = (jj / 2) % 6;

				Oscillator *osc;
				if (jj % 2 == 0)
					osc = new Oscillator(densityHK);
				else
					osc = new Oscillator(densityKD);

				osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

				if (MH == "inverted")
					osc->SetMasses<Oscillator::inverted>(M12, M23);
				else
					osc->SetMasses<Oscillator::normal>(M12, M23);

				TH1D* flux = vFlux.at(jj);
				flux->Add(f1);	//filled with ones

				osc->Oscillate(fin[im], fout[im], flux);

				Reco *reco = copyReco.at(jj);

				for (int ic = 0; ic < 6; ++ic)
				{
					if (chan[ic].find("NC") == std::string::npos)	//it is not a NC
						reco->Scale(chan[ic], flux);
					else if (mode[im] == "nuM0_nuE0" || mode[im] == "nuMB_nuEB") //no NC here
						reco->Scale(chan[ic], 0.0);
					//else no need to scale for NC event
					//reco->Scale(chan[ic], flux);

					TH1D* py = reco->Project(chan[ic], 'y', flux->GetName());

					if (py)
					{
						int bin = nBinAtm + nBinBeam * (offs[ic] + 4*ih);
						if (jj % 2 != 0)
							bin += nBinBeam;

						//std::cout << detName << " " << bin << std::endl;
						for (int n = 0; n < nBinBeam; ++n)
							newmc[bin + n] += py->GetBinContent(n+1);
					}
				}

				//delete flux;
				//delete reco;
				//delete osc;
			}

			for (int jj = 0; jj < copyReco.size(); ++jj)
			{
				delete copyReco.at(jj);
				delete vFlux.at(jj);
			}
	
			newt->Fill();
		}

		std::cout << std::endl;

		newt->Write();
	}


	for (int i = 0; i < vReco.size(); ++i)
		delete vReco.at(i);
	vReco.clear();


	delete f1;

	return 0;
}
