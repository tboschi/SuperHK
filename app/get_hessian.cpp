#include <fstream>
#include <iostream>
#include <chrono>

#include "event/ChiSquared.h"
#include "event/NewChiSquared.h"
#include "physics/Oscillator.h"
#include "physics/ParameterSpace.h"

#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMatrixD.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];
	CardDealer *cd = new CardDealer(cardFile);

	std::string trueOrder, fitOrder;
	cd->Get("true_hierarchy", trueOrder);
	cd->Get("fit_hierarchy", fitOrder);

	//set up oscillation
	Oscillator *osc = new Oscillator(cd);
	//ParameterSpace *parms = new ParameterSpace(cd);

	//create chi2 fitter
	NewChiSquared *fitter = new NewChiSquared(cd);
	fitter->Init();

	int point;
	cd->Get("Point", point);

	//testing at
	//mixing12 0.32 mixing13 0.0256584 mixing23 0.5 deltaCP 0
	//dm12sq 7.6E-5 dm23sq 2.4E-3

	double M12 = 7.6e-5;
	double M23 = 2.4e-3;
	double S12 = 0.32;
	double S13 = 0.0256584;
	double S23 = 0.5;
	double dCP = 0.0;

	//get first file for True point
	//parms->GetEntry(point, M12, M23, S12, S13, S23, dCP);


	//spectra for true point ( On )
	//it is a concat of all events, E_FHC + E_RHC + M_FHC + M_RHC
	if (trueOrder == "normal")
		osc->SetMasses<Oscillator::normal>(M12, M23);
	else if (trueOrder == "inverted")
		osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	Eigen::VectorXd eps = Eigen::VectorXd::Zero(fitter->NumSys());
	Eigen::VectorXd trues = fitter->ConstructSpectrum(osc);
	//Eigen::MatrixXd hessian = fitter->Covariance(trues, trues, eps);
	Eigen::VectorXd gam = fitter->Gamma(trues, eps);

	eps(0) = 1;
	std::vector<double> nchi2 = fitter->ObsX2n(trues, trues, eps);

	std::string output;	//path for out file with extension
	cd->Get("output", output);
	std::ofstream out(output.c_str());
	for (int b = 0; b < nchi2.size(); ++b)
		out << b << "\t" << trues(b) << "\t" << gam(b) 
			 << "\t" << fitter->Scale(gam, eps(0), b)
			 << "\t" << nchi2[b] << std::endl;



	/*
	for (double ss = -0.5; ss <= 0.5; ss += 0.0001) {
		eps(0) = ss;
		out << ss << "\t" << fitter->ObsX2(trues, trues, eps) << "\t"
				  << fitter->SysX2(eps) << std::endl;
	}
	*/


	//std::cout << "\n\n";
	//std::cout << hessian.diagonal() << std::endl;

	//hessian = hessian.inverse();

	//std::string output;	//path for out file with extension
	//cd->Get("output", output);
	//std::ofstream out(output.c_str());
	//for (int x = 0; x < hessian.rows(); ++x)
	//	out << x << "\t" << hessian(x, x) << std::endl;

	out.close();

	delete fitter;

	return 0;
}

