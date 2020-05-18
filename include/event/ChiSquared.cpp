#include "ChiSquared.h"

ChiSquared::ChiSquared(CardDealer *card, int nbins) :
	cd(card),
	_nBin(nbins)
{
	std::string _mode[] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string _chan[] = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	//std::string _type[] = {"E_FHC"};
	std::string _horn[] = {"FHC", "RHC"};
	std::string _oscf[] = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};
	Nu _fin[]  = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	Nu _fout[] = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};

	kMode = sizeof(_mode) / sizeof(std::string);
	kChan = sizeof(_chan) / sizeof(std::string);
	kHorn = sizeof(_horn) / sizeof(std::string);
	kOscf = sizeof(_oscf) / sizeof(std::string);
	kFIn  = sizeof(_fin)  / sizeof(Nu);
	kFOut = sizeof(_fout) / sizeof(Nu);

	mode.insert(mode.end(), &_mode[0], &_mode[kMode]);
	chan.insert(chan.end(), &_chan[0], &_chan[kChan]);
	horn.insert(horn.end(), &_horn[0], &_horn[kHorn]);
	oscf.insert(oscf.end(), &_oscf[0], &_oscf[kOscf]);
	 fin.insert(fin.end(),  &_fin[0],  &_fin[kFIn]  );
	fout.insert(fout.end(), &_fout[0], &_fout[kFOut]);

	if (!cd->Get("sample", type)) {
		std::string _type[] = {"E_FHC", "E_RHC", "M_FHC", "M_RHC"};
		if (nbins < 0)
			kType = sizeof(_type) / sizeof(std::string);
		else
			kType = 1;
		type.insert(type.end(), &_type[0], &_type[kType]);
	}
	else	//sample is defined
		kType = type.size();

	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "types: ";
		for (int t = 0; t < type.size(); ++t)
			std::cout << "\t" << type[t];
		std::cout << std::endl;
	}

	maxIteration = 10;
	err = 1e-9;


	spectrumFile = 0;
	input = "input";
}

void ChiSquared::Init()
{
	LoadCorrelation();
	LoadSystematics();
}

void ChiSquared::LoadCorrelation()
{
	std::string matFile, matName;
	cd->Get("corr_file", matFile);
	cd->Get("corr_name", matName);
	TFile * mF = new TFile(matFile.c_str());
	TMatrixT<double> * cmat = static_cast<TMatrixT<double>*>(mF->Get(matName.c_str()));

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	if (!cd->Get("syst_last", sysB))
		sysB = cmat->GetNcols();

	//corr = Eigen::MatrixXd::Zero(cmat->GetNcols(), cmat->GetNcols());
	corr = Eigen::MatrixXd::Zero(sysB-sysA, sysB-sysA);

	//for (int r = 0; r < cmat->GetNcols(); ++r)
	//	for (int c = 0; c < cmat->GetNcols(); ++c)
	//		corr(r, c) = cmat->operator()(r, c);
	for (int r = sysA; r < sysB; ++r)
		for (int c = sysA; c < sysB; ++c)
			corr(r-sysA, c-sysA) = cmat->operator()(r, c);

	mF->Close();

	corr = corr.inverse();

	if (kVerbosity)
		std::cout << matName << " matrix imported from "
			  << matFile << " with " << sysB-sysA
			  << " entries" << std::endl;
}

/*
 * This routine finds the systematic files and load them into a collection of Matrices.
 * The collection is for spline implementation: it is a mapping between sigma value (int) 
 * and the systematic matrix at that sigma value of the error ( NumBin * NumSys ).
 * If no systematic file is specified in the card, then the matrices are all empty 
 * and it should be like fitting for the stats only, i.e. no fit at all.
 */
void ChiSquared::LoadSystematics()
{
	sysMatrix.clear();
	sysMatrix[0] = Eigen::MatrixXd::Zero(kType * NumBin(), NumSys());
	for (int sigma = -3; sigma < 4; sigma += 2)
		sysMatrix[sigma] = Eigen::MatrixXd::Zero(kType * NumBin(), NumSys());

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	if (!cd->Get("syst_last", sysB))
		sysB = NumSys();

	for (int it = 0; it < kType; ++it)	//type loop
	{
		std::string fileName;
		if (!cd->Get("systematic_" + type[it], fileName))
			continue;

		if (kVerbosity > 1)
			std::cout << "opening file " << fileName << std::endl;
		TFile *sysf = new TFile(fileName.c_str());
		TIter next(sysf->GetListOfKeys());
		TKey *k;
		while (k = static_cast<TKey*>(next()))
		{
			std::string sysname = k->GetName();
			//if (kVerbosity)
				//std::cout << "accepting " << sysname << std::endl;

			int sigma = 0, k_err = 0;
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

				TH1D* hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

				sysname.erase(sysname.find_last_of('_'));
				sysname.erase(0, sysname.find_first_of('_')+1);	//sysname is just number
				k_err = std::stoi(sysname);
				if (k_err < sysA || k_err >= sysB)
					continue;

				//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				for (int n = 0; n < NumBin(); ++n)
					sysMatrix[sigma](it * NumBin() + n, k_err - sysA) = (hsys->GetBinContent(n+1) - 1);
			}
			else
			{
				TH1D* hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

				sysname.erase(0, sysname.find_first_of('_')+1);	//sysname is just number
				k_err = std::stoi(sysname);
				if (k_err < sysA || k_err >= sysB)
					continue;

					//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				for (int n = 0; n < NumBin(); ++n)
					for (sigma = -3; sigma < 4; sigma += 2)
						sysMatrix[sigma](it * NumBin() + n, k_err - sysA) =
								 sigma * (hsys->GetBinContent(n+1) - 1);
			}
		}
		sysf->Close();
	}
}

/*
 * This routine buils the observables as a long vector of NumBin length
 */
Eigen::VectorXd ChiSquared::ConstructSpectrum(Oscillator *osc)
{
	std::string reco_file;

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	std::map<std::string, TH1D*> spectra;
	for (int ih = 0; ih < kHorn; ++ih)	//FHC, RHC
		for (int im = 0; im < kMode; ++im)	//nuE->nuE, nuM->nuE, nuM->nuM
	{
		cd->Get("reco_" + mode[im] + "_" + horn[ih], reco_file);
		Reco *reco = new Reco(reco_file);

		const double *bins;
		int nBin = reco->BinsX(bins);
		f1->SetRange(bins[0], bins[nBin-1]);

		std::string fname = horn[ih] + "_" + mode[im];
		TH1D* flux = new TH1D(fname.c_str(), fname.c_str(), nBin, bins);
		flux->SetDirectory(0);
		flux->Add(f1);	//filled with ones

		osc->Oscillate(fin[im], fout[im], flux);

		for (int ic = 0; ic < kChan; ++ic)	//CCQE, CCnQE, NC x E, M
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
					+ "_" + horn[ih];

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

	double stats;
	if (!cd->Get("stats", stats))
		stats = 1.0;

	Eigen::VectorXd vect(kType * NumBin());
	int j = 0;	//add offset
	for (int it = 0; it < kType; ++it)
	{
		for (int i = 0; i < spectra[type[it]]->GetNbinsX(); ++i, ++j) 
			vect(j) = spectra[type[it]]->GetBinContent(i+1) * stats;

		delete spectra[type[it]];
	}

	return vect;
}

Eigen::VectorXd ChiSquared::LoadSpectrum(int pt, std::string from)
{
	if (kVerbosity)
		std::cout << "passing " << from << std::endl;
	if (spectrumFile)
	{
		if ( (!from.empty() && from != input) ||
		      !PointInFile(spectrumFile, pt) )
		{
			spectrumFile->Close();
			spectrumFile = 0;
		}
	}

	if (!from.empty() && from != input)
		input = from;

	if (!spectrumFile)
	{
		std::string files;
		cd->Get(from, files);

		std::string file;
		std::ifstream listInput(files.c_str());

		while (std::getline(listInput, file))
		{
			spectrumFile = new TFile(file.c_str(), "READ");

			if (spectrumFile->IsZombie())
				continue;

			if (PointInFile(spectrumFile, pt)) //in this file, exit
				break;
			else			  //not in this file, close
			{
				spectrumFile->Close();
				spectrumFile = 0;
			}
		}
	}

	std::cout << "Point from " << spectrumFile->GetName() << std::endl;

	int Point;
	double mcBins[3000];
	TTree * data = static_cast<TTree*>(spectrumFile->Get("mcTree"));
	int bins;

	data->SetBranchAddress("Point",  &Point);
	data->SetBranchAddress("mcBins", mcBins);

	double stats;
	if (!cd->Get("stats", stats))
		stats = 1.0;


	Eigen::VectorXd vect(NumBin());

	for (int i = 0; i < data->GetEntries(); ++ i)
	{
		data->GetEntry(i);

		if (Point == pt)
		{
			for (int n = 0; n < vect.size(); ++n)
				vect(n) = mcBins[n] * stats;
			break;
		}
	}

	return vect;
}

bool ChiSquared::PointInFile(TFile *f, int pt)
{
	int Point;
	double *mcBins;

	TTree * data = static_cast<TTree*>(f->Get("mcTree"));

	data->SetBranchAddress("Point",  &Point);
	data->SetBranchAddress("mcBins", mcBins);

	data->GetEntry(0);
	int pA = Point;
	data->GetEntry(data->GetEntries()-1);
	int pB = Point;

	if (pt >= pA && pt <= pB) //in this file, exit
		return true;
	else
		return false;
}

int ChiSquared::NumSys()
{
	return corr.cols();
}

int ChiSquared::NumBin()
{
	if (_nBin > 0)
		return _nBin;

	std::string reco_file;

	if (!cd->Get("reco_nuE0_nuE0_FHC", reco_file))
		return 0;

	Reco *reco = new Reco(reco_file);
	const double *bins;
	_nBin = reco->BinsY(bins);

	delete reco;

	return _nBin;
}

int ChiSquared::DOF()
{
	return std::abs(kType * NumBin() - NumSys());
}

//On is the true spectrum, En is the observed spectrum
//return time taken for computation
Eigen::VectorXd ChiSquared::FitX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	//initialize epsil with zeroes

	Eigen::VectorXd epsil = Eigen::VectorXd::Zero(NumSys());
	Eigen::VectorXd best_eps = epsil;
	double best_x2 = X2(On, En, best_eps);
	
	int tries = 0;
	while (!FindMinimum(On, En, epsil)/*, alpha)*/
		&& tries < maxIteration)
	{
		++tries;
		double x2 = X2(On, En, epsil);

		if (kVerbosity > 1)
			std::cout << "~~~~~~ No convergence reached ~~~~~~~";

		double dist = (best_eps - epsil).norm();
		double step = x2 - best_x2;
		if (best_x2 > x2)
		{
			best_x2 = x2;
			best_eps = epsil;
			std::cout << "but new best!";
		}


		if (kVerbosity > 1)
		{
			std::cout << " X2 " << best_x2 << std::endl;
			std::cout << "distance from previous best " << dist
				<< " with dx2 " << step << std::endl;
			std::cout << "trying new point " << tries << std::endl;
			std::cout << std::endl;
		}
		epsil.setRandom();
		epsil = best_eps + (std::abs(step) < 1 ? step : 1) * epsil;
	}

	if (tries < maxIteration)
		return epsil;
	else
		return best_eps;
}

//uses Levenberg-Marquardt algorithm
//quites when x2 variation falls below sens and if lambda < 1
//which means it iconvergin
bool ChiSquared::FindMinimum(const Eigen::VectorXd &On, const Eigen::VectorXd &En, 
			     Eigen::VectorXd &epsil)
			     //double alpha)
{
	int c = 0;	//counter
	double lm_0, lm_up, lm_down;	//control parameters
	cd->Get("lm_0", lm_0);
	cd->Get("lm_up", lm_up);
	cd->Get("lm_down", lm_down);

	double lambda = lm_0;

	Eigen::VectorXd delta = Eigen::VectorXd::Constant(NumSys(), 1);

	double x2 = X2(On, En, epsil);
	double diff = 1;
	//double minj = 1;
	//while (diff > err && minj > NumSys() * err && minj > 1e-6 &&
	//       lambda < 1e6 && c < maxIteration)
	while (std::abs(diff) / DOF() > err && delta.norm() / NumSys() > err)
	//while (delta.norm() > NumSys() * err && c < maxIteration)
	{
		++c;

		//build hessian and gradient/jacobian

		//compute first gamma to save computation time
		Eigen::VectorXd gam = Gamma(epsil);

		//corr is inverse of correlation matrix
		Eigen::VectorXd jac = corr * epsil;	//gradient/jacobian
		Eigen::MatrixXd hes = corr;		//hessian

		for (int n = 0; n < kType * NumBin(); ++n)
			for (int k = 0; k < NumSys(); ++k)
			{
				double fkn = Fp(k, n, epsil(k));
				jac(k) += fkn * (En(n) - On(n) / (1 + gam(n)));

				hes(k, k) += On(n) * pow(fkn / (1 + gam(n)), 2);
				for (int j = k+1; j < NumSys(); ++j)
				{
					hes(k, j) += fkn * Fp(j, n, epsil(j))
						   * On(n) / pow(1 + gam(n), 2);
					hes(j, k) = hes(k, j);
				}
			}

		//compute U for stopping criteria
		//add diagonal to hesT hes

		double maxd = hes.diagonal().maxCoeff();
		hes.diagonal() += Eigen::VectorXd::Constant(NumSys(), maxd * lambda);
		delta = hes.ldlt().solve(jac);

		//cos for uphill moves
		//double cosb = delta.dot(newdt) / (newdt.norm() * delta.norm());

		Eigen::VectorXd nextp = epsil - delta;	//next step
		//check if this step is good
		double x2_ = X2(On, En, nextp);
		diff = x2 - x2_;

		//if (x2_ < x2 || (1-cosb) * x2_ < x2)
		if (x2_ < x2)	//next x2 is better, update lambda and epsilons
		{
			//double alpha = 1.0;
			//while (x2 - x2_ < alpha * ctrl * jac.transpose() * delta)
			//{
			//	alpha *= tau;
			//	nextp = epsil - alpha * delta;
			//	x2_ = X2(On, En, nextp);
			//}

			lambda /= lm_down;
			epsil = nextp;
			x2 = x2_;
		}
		else		//next x2 is worse, nothing changes but lambda
			lambda *= lm_up;

		if (kVerbosity > 2)
		//if (false)
		{
			std::cout << c << " -> l " << lambda
				  << ",\tstep: " << std::abs(delta.norm() / NumSys())
				  << ", X2: " << x2
				  << " ( " << std::abs(diff / DOF()) << " ) " << std::endl;
		}
	}


	/*
		//finding best alpha with Armijo-Goldstein condition
		double alpha = 1.0;
		Eigen::VectorXd nextp = epsil - alpha * delta;
		while (X2(On, En, epsil) - X2(On, En, nextp) 
			< alpha * ctrl * jac.transpose() * delta)
		{
			alpha *= tau;
			nextp = epsil - alpha * delta;
		}
	}
	*/

	if (lambda < 1./lm_0)	//convergence was reached
		return true;
	else			//no convergence, change starting point
		return false;
}

Eigen::MatrixXd ChiSquared::Covariance(const Eigen::VectorXd &On, const Eigen::VectorXd &En, const Eigen::VectorXd &epsil)
{
	// scaling by X2 / dof
	//double res = X2(On, En, epsil) / DOF();

	//compute hessian
	Eigen::VectorXd gam = Gamma(epsil);
	Eigen::MatrixXd hes = corr;
	for (int n = 0; n < kType * NumBin(); ++n) {
		for (int k = 0; k < NumSys(); ++k) {
			double fkn = Fp(k, n, epsil(k));
			hes(k, k) += On(n) * pow(fkn / (1 + gam(n)), 2);
			for (int j = k + 1; j < NumSys(); ++j) {
				hes(k, j) += fkn * Fp(j, n, epsil(j))
					   * On(n) / pow(1 + gam(n), 2);
				hes(j, k) = hes(k, j);
			}
		}
	}

	return hes.inverse();
	//return res * hes.inverse();
}


std::string ChiSquared::Diag(const Eigen::VectorXd &On)
{

	//return res * hes.inverse();
	Eigen::MatrixXd hes = Eigen::MatrixXd::Zero(NumSys(), NumSys());

	//std::cout << "corr\n" << hes << std::endl;
	for (double s = -2; s < 2.1; s += 0.2) {
		double sum = 0;
		for (int n = 0; n < kType * NumBin(); ++n)
			sum += On(n) * (1 + F(9, n, s));
		std::cout << s << "\t" << sum << std::endl;
	}

	//std::cout << "corr\n" << hes << std::endl;
	for (int n = 0; n < kType * NumBin(); ++n)
		for (int k = 0; k < NumSys(); ++k)
		{
			double fkn = Fp(k, n, 0);
			hes(k, k) += On(n) * fkn * fkn;
			for (int j = k + 1; j < NumSys(); ++j) {
				hes(k, j) += On(n) * fkn * Fp(j, n, 0);
				hes(j, k) = hes(k, j);
			}
		}

	std::stringstream ss;
	for (int n = 0; n < kType * NumBin(); ++n) {
		ss << On(n);
		for (int k = 0; k < NumSys(); ++k)
			ss << "\t" << Fp(k, n, 0);
		ss << "\n";
	}

	//ss << hes.diagonal() << "\n\n";
	//hes += corr;
	//ss << hes.inverse().diagonal() << std::endl;

	return ss.str();
}


/*
//tries to define a contour for the kth systematic
void ChiSquared::DefineContour(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			       double alpha, double &e_min, double &e_max,
			       int k_err, double x2min)
{
	//make a copy of epsil
	Eigen::VectorXd backup = epsil;

	int c = 0;
	double err = 1e-9;

	std::vector<double> eLimits(2);
	//two iterations, left and right of e_k
	for (int j = 0; j < 2; ++j)
	{
		//shift e_k by 10%, multiplying by 0.9 (j == 0) or 1.1 (j == 1)
		epsil(k_err) *= (0.9 + j * 0.2);

		Eigen::VectorXd delta = Eigen::VectorXd::Constant(NumSys(), 1);
		while (delta.norm() > err && c < maxIterations)
		{
			++c;

			//compute first gamma to save computation time
			Eigen::VectorXd gam = Gamma();

			//corr is inverse of correlation matrix
			Eigen::VectorXd der = corr * epsil;
			Eigen::MatrixXd jac = corr;

			for (int n = 0; n < kType * NumBin(); ++n)
			{
				for (int k = 0; k < NumSys(); ++k)
				{
					der(k) += Fp(k, n) * (En(n) - On(n) / (1 + gam(n)));

					for (int j = 0; j < NumSys(); ++j)
						jac(j, k) += Fp(k, n) * Fp(j, n) * On(n) /
							     pow(1 + gam(n), 2);
				}
			}

			jac.row(k_err) = der.transpose();
			der(k_err) = X2(On, En) - x2min - 1;

			//using Norm eq.
			delta = (jac.transpose() * jac).ldlt().solve(jac.transpose() * der);

			epsil -= alpha * delta;		//alpha regulates convergence
		}

		eLimits[j] = epsil(k_err);

		//restore epsil
		epsil = backup;
	}

	if (eLimits[0] <= eLimits[1])
	{
		e_min = eLimits[0];
		e_max = eLimits[1];
	}
	else
	{
		e_min = eLimits[1];
		e_max = eLimits[0];
	}
}
*/

double ChiSquared::X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		      const Eigen::VectorXd &epsil)
{
	return ObsX2(On, En, epsil) + SysX2(epsil);
}

/*
double ChiSquared::X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	return X2(On, En, epsil);
}
*/

/*
double ChiSquared::ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	return ObsX2(On, En, epsil);
}
*/

double ChiSquared::ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			 const Eigen::VectorXd &epsil)
{
	//if no systematic, gamma is just zeros
	Eigen::VectorXd gam = Gamma(epsil);

	double chi2 = 0;
	for (int n = 0; n < kType * NumBin(); ++n)
	{
		//std::cout << "obs " << 1+gam(n) << "\t" << En(n) << "\t" << On(n)
		//	  << "\t" << On(n) * log((1+gam(n)) * En(n) / On(n)) << std::endl;
		chi2 += 2 * ((1 + gam(n)) * En(n) - On(n));
		if (1 + gam(n) > 1e-9 && On(n) > 1e-9 && En(n) > 1e-9)
			chi2 -= 2 * On(n) * log((1 + gam(n)) * En(n) / On(n));
	}

	return chi2;
}

double ChiSquared::SysX2(const Eigen::VectorXd &epsil)
{
	return epsil.transpose() * corr * epsil;
}

/*
Eigen::VectorXd ChiSquared::Gamma()
{
	return Gamma(epsil);
}
*/

Eigen::VectorXd ChiSquared::Gamma(const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd gam = Eigen::VectorXd::Zero(kType * NumBin());
	for (int n = 0; n < kType * NumBin(); ++n)
		for (int k = 0; k < NumSys(); ++k)
			gam(n) += F(k, n, epsil(k));

	return gam;
}

/*
double ChiSquared::F(int k, int n)
{
	return F(k, n, epsil);
}
*/

double ChiSquared::F(int k, int n, double eij)
{
	double dl, du;

	if (eij < -1)
	{
		dl = -3; du = -1;
	}
	else if (eij >= -1 && eij < 0)
	{
		dl = -1; du = 0;
	}
	else if (eij >= 0 && eij < 1)
	{
		dl = 0; du = 1;
	}
	else if (eij >= 1)
	{
		dl = 1; du = 3;
	}

	return Fp(k, n, dl, du) * (eij - dl) + sysMatrix[dl](n, k);
	//return F(k, n, dl, du, eij);
}

/*
double ChiSquared::F(int k, int n, double dl, double du)
{
	return F(k, n, dl, du, epsil);
}
*/
/*
double ChiSquared::F(int k, int n, double dl, double du, double eij)
{
	return Fp(k, n, dl, du) * (eij - dl) + sysMatrix[dl](n, k);
}
*/

/*
double ChiSquared::Fp(int k, int n)
{
	return Fp(k, n, epsil);
}
*/

double ChiSquared::Fp(int k, int n, double eij)
{
	double dl, du;

	if (eij < -1)
	{
		dl = -3; du = -1;
	}
	else if (eij >= -1 && eij < 0)
	{
		dl = -1; du = 0;
	}
	else if (eij >= 0 && eij < 1)
	{
		dl = 0; du = 1;
	}
	else if (eij >= 1)
	{
		dl = 1; du = 3;
	}

	return Fp(k, n, dl, du);
}

double ChiSquared::Fp(int k, int n, double dl, double du)
{
	return (sysMatrix[du](n, k) - sysMatrix[dl](n, k)) / (du - dl);
}

/*
Eigen::VectorXd ChiSquared::Epsilons()
{
	return epsil;
}
*/
