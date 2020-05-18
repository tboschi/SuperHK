#include <fstream>
#include <iostream>
#include <iomanip>

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

	std::string treeName, files, outname;
	double dvf;
	cd->Get("tree", treeName);
	cd->Get("files", files);
	cd->Get("output", outname);
	cd->Get("division", dvf);

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

	TChain *ch = new TChain(treeName.c_str());
	ch->Add(files.c_str());

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

	Oscillator *oscHK = new Oscillator(densityHK);
	Oscillator *oscKD = new Oscillator(densityKD);

	ch->SetBranchAddress("mcBins", mcBins);
	ch->SetBranchAddress("Point", &Point);
	ch->SetBranchAddress("CP",  &dCP); 
	ch->SetBranchAddress("M12", &M12);
	ch->SetBranchAddress("M23", &M23);
	ch->SetBranchAddress("S12", &S12);
	ch->SetBranchAddress("S13", &S13);
	ch->SetBranchAddress("S23", &S23);

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);
	f1->SetRange(bins[0], bins[nBinBeam-1]);

	std::stringstream ssl;
	ssl << nBinAtm + nBinBeam * 8;
	std::string mcbin = "mcBins[" + ssl.str() + "]/D";
	//std::cout << "allocating " << ssl.str() << std::endl;

	//TObjArray *fileElements = ch->GetListOfFiles();
	//TIter next(fileElements);
	//TChainElement *current = 0;
	//int nentry = 0;
	//ch->GetEntry(0);	//load first entry to compute nentries
	//while ( current = static_cast<TChainElement*>(next()) )
	//{
	//	std::string name(current->GetTitle());
	//	//std::cout << "Reading " << name << std::endl;
	//	name.insert(name.find(".root"), "_oscillated");
	//	TFile * outf = new TFile(name.c_str(), "RECREATE");
	//	outf->cd();
	//}

	int ntot = ch->GetEntries();
	int cc = ntot / int(dvf);
	int rr = ntot % int(dvf);

	int nentry = 0;
	int ww = log10(dvf) + 1;
	for (int f = 0; f < dvf; ++f)
	{
		//std::stringstream ssn;
		//ssn << "_" << std::setfill('0') << std::setw(ww) << f;
		std::string name = outname;
		name.insert(name.find(".000"), "_oscillated");
		std::cout << "creating file " << name << std::endl;
		TFile *outf = new TFile(name.c_str(), "RECREATE");
		outf->cd();

		TTree *newt = new TTree("mcTree", ch->GetTitle());

		newt->Branch("mcBins", newmc, mcbin.c_str());
		newt->Branch("Point", &Point, "Point/I");
		newt->Branch("CP",  &dCP, "CP/D"); 
		newt->Branch("M12", &M12, "M12/D");
		newt->Branch("M23", &M23, "M23/D");
		newt->Branch("S12", &S12, "S12/D");
		newt->Branch("S13", &S13, "S13/D");
		newt->Branch("S23", &S23, "S23/D");

		double _dCP = -9999;
		int evt = (f < rr) ? (cc + 1) : cc;
		for (int n = nentry; n < nentry + evt && n < ntot; ++n)
		{
			ch->GetEntry(n);

			if (dCP != _dCP)
			{
				std::cout << "dCP --- " << dCP << std::endl;
				_dCP = dCP;
			}

			//copy atmospheric
			memcpy(newmc, mcBins,      sizeof(double) * nBinAtm);
			memset(newmc + nBinAtm, 0, sizeof(double) * nBinBeam * 8);

			//double number[4];
			//for (int i = 0; i < 4; ++i)
			//	number[i] = 0;
			//generate beam

			oscHK->Reset();
			oscKD->Reset();

			oscHK->SetMasses<Oscillator::normal>(M12, M23);
			oscKD->SetMasses<Oscillator::normal>(M12, M23);

			oscHK->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
			oscKD->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

			for (int ih = 0; ih < 2; ++ih)
			{
				for (int im = 0; im < 6; ++im)
				{
					int jj = im + 6 * ih;

					std::string fnameHK = horn[ih] + "_" + mode[im] + "_HK";
					std::string fnameKD = horn[ih] + "_" + mode[im] + "_KD";

					TH1D* fluxHK = new TH1D(fnameHK.c_str(), fnameHK.c_str(), nBinEnu, bins);
					TH1D* fluxKD = new TH1D(fnameKD.c_str(), fnameKD.c_str(), nBinEnu, bins);

					fluxHK->SetDirectory(0);
					fluxKD->SetDirectory(0);

					fluxHK->Add(f1);	//filled with ones
					fluxKD->Add(f1);	//filled with ones

					oscHK->Oscillate(fin[im], fout[im], fluxHK);
					oscKD->Oscillate(fin[im], fout[im], fluxKD);

					Reco *recoHK = new Reco(*vReco.at(jj));
					Reco *recoKD = new Reco(*vReco.at(jj));

					for (int ic = 0; ic < 6; ++ic)
					{
						if (chan[ic].find("NC") == std::string::npos)	//it is not a NC
						{
							recoHK->Scale(chan[ic], fluxHK);
							recoKD->Scale(chan[ic], fluxKD);
						}
						else if (mode[im] == "nuM0_nuE0" || mode[im] == "nuMB_nuEB") //no NC here
						{
							recoHK->Scale(chan[ic], 0.0);
							recoKD->Scale(chan[ic], 0.0);
						}
						//else no need to scale for NC event
						//reco->Scale(chan[ic], flux);

						TH1D* pyHK = recoHK->Project(chan[ic], 'y', "_HK");
						TH1D* pyKD = recoKD->Project(chan[ic], 'y', "_KD");

						if (pyHK && pyKD)
						{
							int binHK = nBinAtm + nBinBeam * (offs[ic] + 4*ih + 0);
							int binKD = nBinAtm + nBinBeam * (offs[ic] + 4*ih + 1);

							for (int n = 0; n < nBinBeam; ++n)
							{
								newmc[binHK + n] += pyHK->GetBinContent(n+1);
								newmc[binKD + n] += pyKD->GetBinContent(n+1);
							}
						}
					}

					delete fluxHK;
					delete fluxKD;
					delete recoHK;
					delete recoKD;
				}
			}
			newt->Fill();
		}

		nentry += evt;

		newt->Write();
		outf->Close();
	}

	for (int i = 0; i < vReco.size(); ++i)
		delete vReco.at(i);
	vReco.clear();


	delete f1;
	delete ch;

	return 0;
}
