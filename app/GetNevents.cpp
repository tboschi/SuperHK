#include <iostream>
#include <iomanip>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
//#include "TIter.h"
#include "TKey.h"
#include "TClass.h"

int main(int argc, char** argv)
{
	TFile* inF = new TFile(argv[1], "OPEN");
	//TFile* inR = new TFile(argv[2], "OPEN");


	std::map<std::string, double> mInt;
	std::map<std::string, double>::iterator im;

	TIter next(inF->GetListOfKeys());
	TKey *key;
	std::cout << "Keys found in " << argv[1] << std::endl;
	while ( key = static_cast<TKey*>(next()) )
	{
		if (key->GetClassName() == "TH1D");
		{
			TH1D *h = static_cast<TH1D*>(key->ReadObj());
			std::cout << "\t> " << h->GetName() << std::endl;
			mInt[std::string(h->GetName())] = h->Integral("");
		}
	}

	/*		ORIGINAL
	double nuE0_nuE0_E_CCQE  = static_cast<TH1D*>(inF->Get("nuE0_nuE0_E_CCQE"))->Integral("");
	double nuE0_nuE0_E_CCnQE = static_cast<TH1D*>(inF->Get("nuE0_nuE0_E_CCnQE"))->Integral("");
	double nuE0_nuE0_E_NC    = static_cast<TH1D*>(inF->Get("nuE0_nuE0_E_NC"))->Integral("");
	double nuE0_nuE0_M_CCQE  = static_cast<TH1D*>(inF->Get("nuE0_nuE0_M_CCQE"))->Integral("");
	double nuE0_nuE0_M_CCnQE = static_cast<TH1D*>(inF->Get("nuE0_nuE0_M_CCnQE"))->Integral("");
	double nuE0_nuE0_M_NC    = static_cast<TH1D*>(inF->Get("nuE0_nuE0_M_NC"))->Integral("");
	double nuM0_nuM0_E_CCQE  = static_cast<TH1D*>(inF->Get("nuM0_nuM0_E_CCQE"))->Integral("");
	double nuM0_nuM0_E_CCnQE = static_cast<TH1D*>(inF->Get("nuM0_nuM0_E_CCnQE"))->Integral("");
	double nuM0_nuM0_E_NC    = static_cast<TH1D*>(inF->Get("nuM0_nuM0_E_NC"))->Integral("");
	double nuM0_nuM0_M_CCQE  = static_cast<TH1D*>(inF->Get("nuM0_nuM0_M_CCQE"))->Integral("");
	double nuM0_nuM0_M_CCnQE = static_cast<TH1D*>(inF->Get("nuM0_nuM0_M_CCnQE"))->Integral("");
	double nuM0_nuM0_M_NC    = static_cast<TH1D*>(inF->Get("nuM0_nuM0_M_NC"))->Integral("");
	double nuM0_nuE0_E_CCQE  = static_cast<TH1D*>(inF->Get("nuM0_nuE0_E_CCQE"))->Integral("");
	double nuM0_nuE0_E_CCnQE = static_cast<TH1D*>(inF->Get("nuM0_nuE0_E_CCnQE"))->Integral("");
	double nuM0_nuE0_E_NC    = static_cast<TH1D*>(inF->Get("nuM0_nuE0_E_NC"))->Integral("");
	double nuM0_nuE0_M_CCQE  = static_cast<TH1D*>(inF->Get("nuM0_nuE0_M_CCQE"))->Integral("");
	double nuM0_nuE0_M_CCnQE = static_cast<TH1D*>(inF->Get("nuM0_nuE0_M_CCnQE"))->Integral("");
	double nuM0_nuE0_M_NC    = static_cast<TH1D*>(inF->Get("nuM0_nuE0_M_NC"))->Integral("");

	double nuEB_nuEB_E_CCQE  = static_cast<TH1D*>(inF->Get("nuEB_nuEB_E_CCQE"))->Integral("");
	double nuEB_nuEB_E_CCnQE = static_cast<TH1D*>(inF->Get("nuEB_nuEB_E_CCnQE"))->Integral("");
	double nuEB_nuEB_E_NC    = static_cast<TH1D*>(inF->Get("nuEB_nuEB_E_NC"))->Integral("");
	double nuEB_nuEB_M_CCQE  = static_cast<TH1D*>(inF->Get("nuEB_nuEB_M_CCQE"))->Integral("");
	double nuEB_nuEB_M_CCnQE = static_cast<TH1D*>(inF->Get("nuEB_nuEB_M_CCnQE"))->Integral("");
	double nuEB_nuEB_M_NC    = static_cast<TH1D*>(inF->Get("nuEB_nuEB_M_NC"))->Integral("");
	double nuMB_nuMB_E_CCQE  = static_cast<TH1D*>(inF->Get("nuMB_nuMB_E_CCQE"))->Integral("");
	double nuMB_nuMB_E_CCnQE = static_cast<TH1D*>(inF->Get("nuMB_nuMB_E_CCnQE"))->Integral("");
	double nuMB_nuMB_E_NC    = static_cast<TH1D*>(inF->Get("nuMB_nuMB_E_NC"))->Integral("");
	double nuMB_nuMB_M_CCQE  = static_cast<TH1D*>(inF->Get("nuMB_nuMB_M_CCQE"))->Integral("");
	double nuMB_nuMB_M_CCnQE = static_cast<TH1D*>(inF->Get("nuMB_nuMB_M_CCnQE"))->Integral("");
	double nuMB_nuMB_M_NC    = static_cast<TH1D*>(inF->Get("nuMB_nuMB_M_NC"))->Integral("");
	double nuMB_nuEB_E_CCQE  = static_cast<TH1D*>(inF->Get("nuMB_nuEB_E_CCQE"))->Integral("");
	double nuMB_nuEB_E_CCnQE = static_cast<TH1D*>(inF->Get("nuMB_nuEB_E_CCnQE"))->Integral("");
	double nuMB_nuEB_E_NC    = static_cast<TH1D*>(inF->Get("nuMB_nuEB_E_NC"))->Integral("");
	double nuMB_nuEB_M_CCQE  = static_cast<TH1D*>(inF->Get("nuMB_nuEB_M_CCQE"))->Integral("");
	double nuMB_nuEB_M_CCnQE = static_cast<TH1D*>(inF->Get("nuMB_nuEB_M_CCnQE"))->Integral("");
	double nuMB_nuEB_M_NC    = static_cast<TH1D*>(inF->Get("nuMB_nuEB_M_NC"))->Integral("");
	*/

	//double RnuE0_nuE0_E_CCQE  = static_cast<TH1D*>(inR->Get("nuE0_nuE0_E_CCQE"))->Integral("");
	//double RnuE0_nuE0_E_CCnQE = static_cast<TH1D*>(inR->Get("nuE0_nuE0_E_CCnQE"))->Integral("WIDTH");
	//double RnuE0_nuE0_E_NC    = static_cast<TH1D*>(inR->Get("nuE0_nuE0_E_NC"))->Integral("WIDTH");
	//double RnuE0_nuE0_M_CCQE  = static_cast<TH1D*>(inR->Get("nuE0_nuE0_M_CCQE"))->Integral("WIDTH");
	//double RnuE0_nuE0_M_CCnQE = static_cast<TH1D*>(inR->Get("nuE0_nuE0_M_CCnQE"))->Integral("WIDTH");
	//double RnuE0_nuE0_M_NC    = static_cast<TH1D*>(inR->Get("nuE0_nuE0_M_NC"))->Integral("WIDTH");
	//double RnuM0_nuM0_E_CCQE  = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_CCQE"))->Integral("WIDTH");
	//double RnuM0_nuM0_E_CCnQE = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_CCnQE"))->Integral("WIDTH");
	//double RnuM0_nuM0_E_NC    = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_NC"))->Integral("WIDTH");
	//double RnuM0_nuM0_M_CCQE  = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_CCQE"))->Integral("WIDTH");
	//double RnuM0_nuM0_M_CCnQE = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_CCnQE"))->Integral("WIDTH");
	//double RnuM0_nuM0_M_NC    = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_NC"))->Integral("WIDTH");
	//double RnuM0_nuE0_E_CCQE  = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_CCQE"))->Integral("WIDTH");
	//double RnuM0_nuE0_E_CCnQE = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_CCnQE"))->Integral("WIDTH");
	//double RnuM0_nuE0_E_NC    = static_cast<TH1D*>(inR->Get("nuM0_nuE0_E_NC"))->Integral("WIDTH");
	//double RnuM0_nuE0_M_CCQE  = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_CCQE"))->Integral("WIDTH");
	//double RnuM0_nuE0_M_CCnQE = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_CCnQE"))->Integral("WIDTH");
	//double RnuM0_nuE0_M_NC    = static_cast<TH1D*>(inR->Get("nuM0_nuE0_M_NC"))->Integral("WIDTH");

	//double RnuEB_nuEB_E_CCQE  = static_cast<TH1D*>(inR->Get("nuEB_nuEB_E_CCQE"))->Integral("WIDTH");
	//double RnuEB_nuEB_E_CCnQE = static_cast<TH1D*>(inR->Get("nuEB_nuEB_E_CCnQE"))->Integral("WIDTH");
	//double RnuEB_nuEB_E_NC    = static_cast<TH1D*>(inR->Get("nuEB_nuEB_E_NC"))->Integral("WIDTH");
	//double RnuEB_nuEB_M_CCQE  = static_cast<TH1D*>(inR->Get("nuEB_nuEB_M_CCQE"))->Integral("WIDTH");
	//double RnuEB_nuEB_M_CCnQE = static_cast<TH1D*>(inR->Get("nuEB_nuEB_M_CCnQE"))->Integral("WIDTH");
	//double RnuEB_nuEB_M_NC    = static_cast<TH1D*>(inR->Get("nuEB_nuEB_M_NC"))->Integral("WIDTH");
	//double RnuMB_nuMB_E_CCQE  = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_CCQE"))->Integral("WIDTH");
	//double RnuMB_nuMB_E_CCnQE = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_CCnQE"))->Integral("WIDTH");
	//double RnuMB_nuMB_E_NC    = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_NC"))->Integral("WIDTH");
	//double RnuMB_nuMB_M_CCQE  = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_CCQE"))->Integral("WIDTH");
	//double RnuMB_nuMB_M_CCnQE = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_CCnQE"))->Integral("WIDTH");
	//double RnuMB_nuMB_M_NC    = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_NC"))->Integral("WIDTH");
	//double RnuMB_nuEB_E_CCQE  = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_CCQE"))->Integral("WIDTH");
	//double RnuMB_nuEB_E_CCnQE = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_CCnQE"))->Integral("WIDTH");
	//double RnuMB_nuEB_E_NC    = static_cast<TH1D*>(inR->Get("nuMB_nuEB_E_NC"))->Integral("WIDTH");
	//double RnuMB_nuEB_M_CCQE  = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_CCQE"))->Integral("WIDTH");
	//double RnuMB_nuEB_M_CCnQE = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_CCnQE"))->Integral("WIDTH");
	//double RnuMB_nuEB_M_NC    = static_cast<TH1D*>(inR->Get("nuMB_nuEB_M_NC"))->Integral("WIDTH");

	std::ofstream tout (argv[2]);
	tout << std::setfill(' ') << std::setprecision(5) << std::fixed;
	tout << "\\textbf{1 Ring $\\nu_\\mu$-like}\n";
	//1 Ring nu_mu like
	tout << "\\begin{tabular}{lrrrrrrr}" << "\n";
	tout << "\t\\toprule" << "\n";
	tout << "\t& $\\nu_\\mu$"
	     << "\t& $\\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu$"
	     << "\t& $\\cj{\\nu}_e$"
	     << "\t& $\\nu_\\mu \\to \\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu \\to \\cj{\\nu}_e$"
	     << "\t& Total"
	     << "\t\\\\" << "\n";
	tout << "\t\t\\midrule" << "\n";
	//CCQE numbers
	tout << "\tCCQE\t& "
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
	tout << "\tCCnQE\t& "
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
	tout << "\tNC\t& "
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
	tout << "\t\\midrule" << "\n";	
	//Total numbers
	tout << "\tTot\t& "
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
	tout << "\t\\bottomrule" << "\n";	
	tout << "\\end{tabular}" << "\n";

	tout << "\n\\vfill\n";

	tout << "\\textbf{1 Ring $\\nu_e$-like}\n";
	//1 Ring nu_e like
	tout << "\\begin{tabular}{lrrrrrrr}" << "\n";
	tout << "\t\\toprule" << "\n";
	tout << "\t& $\\nu_\\mu$"
	     << "\t& $\\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu$"
	     << "\t& $\\cj{\\nu}_e$"
	     << "\t& $\\nu_\\mu \\to \\nu_e$"
	     << "\t& $\\cj{\\nu}_\\mu \\to \\cj{\\nu}_e$"
	     << "\t& Total"
	     << "\t\\\\" << "\n";
	tout << "\t\\midrule" << "\n";
	//CCQE numbers
	tout << "\tCCQE\t& "
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
	tout << "\tCCnQE\t& "
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
	tout << "\tNC\t& "
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
	tout << "\t\\midrule" << "\n";	
	//Total numbers
	tout << "\tTot\t& "
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
	tout << "\t\\bottomrule" << "\n";	
	tout << "\\end{tabular}" << "\n";

	return 0;
}
