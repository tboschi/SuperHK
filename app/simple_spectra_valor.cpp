#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"

#include "valor/Conventions/Constants.h"
#include "Binning.h"

using namespace valor;
using namespace valor::util;
using namespace valor::conv;
using namespace std;

void simple_spectra(bool oscinput = false, const char * outtag="_valorosc", bool noosc=false, 
		const char * osc3ppinputloc="."
		//const char * osc3ppinputloc="/home/hyperk/tdealtry/public_html/oa/190404/"
		)
{
	string str_samples[2] = {"1Rmu","1Re"};
	string str_modes[3] = {"CCQE","CCnQE","NC"};
	string infiletags[6] = {"numu_x_numu", "nue_x_nue",
		"numubar_x_numubar", "nuebar_x_nuebar",
		"numu_x_nue", "numubar_x_nuebar"};
	string infiletags_rhc[6] = {"numu_crs_numu_anu", "nue_crs_nue_anu",
		"numubar_crs_numubar_anu", "nuebar_crs_nuebar_anu",
		"numu_crs_nue_anu", "numubar_crs_nuebar_anu"};  
	string oschistname[6] = {"oscMM_0", "oscEE_0",
		"oscMM_B", "oscEE_B",
		"oscME_0", "oscME_B"};
	TFile * fevt = 0;
	TFile fosc(oscinput ? "osctest.root" : "valoroscprob.root");
	TH1D * hosc = 0;
	TH2D * hevt = 0;
	ofstream tf;
	tf.open(TString::Format("osctable%s.tex", outtag));
	tf << fixed << setprecision(5);
	tf << "\\documentclass[10pt,a4paper]{article}" << endl
		<< "\\usepackage[utf8]{inputenc}" << endl
		<< "\\usepackage{amsmath}" << endl
		<< "\\usepackage{amsfonts}" << endl
		<< "\\usepackage{amssymb}" << endl
		<< "\\usepackage{fullpage}" << endl
		<< "\\begin{document}" << endl
		<< "\\begin{center}" << endl;
	TFile fout("simple.root", "RECREATE");
	for(int is = 0; is < 4; is++) {
		SuperkSample_t s = kSKSamples[is];
		double n_in_sample = 0;
		tf << "\\begin{tabular}{|l|r|r|r|r|r|r||r|}" << endl
			<< "\\hline" << endl
			<< (is < 2 ? "FHC" : "RHC") << str_samples[is%2]
			<< " & $\\nu_\\mu$ & $\\nu_e$ & $\\bar{\\nu}_\\mu$ & $\\bar{\\nu}_e$ & Osc. $\\nu_e$ & Osc. $\\bar{\\nu}_e$ & Total \\\\" << endl
			<< "\\hline" << endl;
		TH2D *hnewsample = new TH2D(TString::Format("enu_erec_%s%s_TOTAL", (is < 2 ? "" : "RHC"), str_samples[is%2].c_str()),
				";Etrue;Ereco",
				kNEtrueBins[is], util::EtrueBinEdges(s),
				kNErecoBins[is], util::ErecoBinEdges(s));
		for(int im = 0; im < 3; im++) {
			double n_in_mode = 0;
			tf << str_modes[im];
			TH2D * hnew = new TH2D(TString::Format("enu_erec_%s%s_%s", (is < 2 ? "" : "RHC"), str_samples[is%2].c_str(), str_modes[im].c_str()),
					";Etrue;Ereco",
					kNEtrueBins[is], util::EtrueBinEdges(s),
					kNErecoBins[is], util::ErecoBinEdges(s));
			for(int io = 0; io < 6; io++) {
				if((im == 2 || noosc) && io == 4) {
					tf << " & &";
					break; //skip oscillated NC 
				}
				double n_in_oscmode = 0;
				fevt = new TFile(TString::Format("%s/enu_erec_%s.root", osc3ppinputloc, is < 2 ? infiletags[io].c_str() : infiletags_rhc[io].c_str()));
				fosc.GetObject(oschistname[io].c_str(), hosc);
				fevt->GetObject(TString::Format("enu_erec_%s_%s", str_samples[is%2].c_str(), str_modes[im].c_str()), hevt);
				for(int ireco = 1; ireco <= hevt->GetNbinsY(); ireco++) {
					double n_in_recobin = 0;
					for(int itrue = 1; itrue <= hevt->GetNbinsX(); itrue++) {
						double oscprob = ((im == 2 || noosc) ? 1 : hosc->GetBinContent(itrue)); //don't apply oscillation to NC events
						double rawevts = hevt->GetBinContent(itrue, ireco);
						double oscevts = oscprob * rawevts;
						n_in_recobin += oscevts;
						hnew->Fill(util::EtrueBinEdges(s,itrue-1)+conv::kOneMeV,
								util::ErecoBinEdges(s,ireco-1)+conv::kOneMeV,
								oscevts);
						hnewsample->Fill(util::EtrueBinEdges(s,itrue-1)+conv::kOneMeV,
								util::ErecoBinEdges(s,ireco-1)+conv::kOneMeV,
								oscevts);
					}//itrue
					//cout << "RECOBIN " << hevt->GetYaxis()->GetBinCenter(ireco) << " \t " << n_in_recobin << endl;
					n_in_oscmode += n_in_recobin;
				}//ireco
				fevt->Close();
				cout << "OSCMODE " << (is < 2 ? infiletags[io] : infiletags_rhc[io].c_str()) << " \t " << n_in_oscmode << endl;
				n_in_mode += n_in_oscmode;
				tf << " & " << n_in_oscmode;
			}//io
			cout << "MODE " << str_modes[im] << " \t " << n_in_mode << endl;
			n_in_sample += n_in_mode;
			tf << " & " << n_in_mode << " \\\\" << endl;
			fout.cd();
			hnew->Write();
			delete hnew;
		}//im
		cout << "SAMPLE " << (is < 2 ? "FHC" : "RHC") << str_samples[is%2] << " \t " << n_in_sample << endl << endl;
		tf << "\\hline" << endl
			<< "Total & &&&&&& " << n_in_sample << " \\\\" << endl
			<< "\\hline" << endl
			<< "\\end{tabular}" << endl << endl
			<< "\\vspace{2cm}" << endl << endl;
		hnewsample->Write();
		delete hnewsample;
	}//is
	tf << "\\end{center}" << endl
		<< "\\end{document}" << endl;
	tf.close();
}

