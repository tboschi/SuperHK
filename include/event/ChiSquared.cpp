#include "NewChiSquared.h"

NewChiSquared::NewChiSquared(CardDealer *card, int nbins) :
	cd(card),
	_nBin(nbins)
{
	std::string _mode[] = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0",
			       "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	std::string _chan[] = {"E_CCQE", "M_CCQE", "E_CCnQE",
			       "M_CCnQE", "E_NC", "M_NC"};
	std::string _horn[] = {"FHC", "RHC"};
	std::string _oscf[] = {"oscEE_0", "oscMM_0", "oscME_0",
			       "oscEE_B", "oscMM_B", "oscME_B"};
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

	if (!cd->Get("scale_", scale)) {
		double err;
		scale.clear();
		if (!cd->Get("scale", err))
			scale[""] = 0.024;
		else
			scale[""] = err;
	}
	kScale = scale.size();	//of size 1, 2, or 4


	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "types: ";
		for (int t = 0; t < type.size(); ++t)
			std::cout << "\t" << type[t];
		std::cout << std::endl;
		std::cout << "scales: " << kScale << std::endl;
		for (is = scale.begin(); is != scale.end(); ++is) 
			std::cout << "\t" << is->first << " -> " << is->second << std::endl;
	}

	maxIteration = 10;
	err = 1e-9;

	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	spectrumFile = 0;
	input = "input";
}

void NewChiSquared::Init()
{
	LoadCorrelation();
	LoadSystematics();
}

void NewChiSquared::LoadCorrelation()
{
	std::string matFile, matName;
	cd->Get("corr_file", matFile);
	cd->Get("corr_name", matName);
	TFile * mf = new TFile(matFile.c_str());
	TMatrixT<double> * cmat = static_cast<TMatrixT<double>*>(mf->Get(matName.c_str()));

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	if (!cd->Get("syst_last", sysB))
		sysB = cmat->GetNcols();

	if (kVerbosity)
		std::cout << "Limiting systematics between "
			  << sysA << " and " << sysB << std::endl;
	//first kScale entries are for energy scale
	corr = Eigen::MatrixXd::Zero(kScale + sysB - sysA, kScale + sysB - sysA);
	corr.topLeftCorner(kScale, kScale) =
				Eigen::MatrixXd::Identity(kScale, kScale);

	for (int r = sysA; r < sysB; ++r)
		for (int c = sysA; c < sysB; ++c)
			corr(r + kScale - sysA, c + kScale - sysA) = (*cmat)(r, c);
	mf->Close();

	corr = corr.inverse();
	std::cout << "corr diag\n" << corr.diagonal() << std::endl;


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
void NewChiSquared::LoadSystematics()
{
	sysMatrix.clear();
	sysMatrix[0] = Eigen::MatrixXd::Zero(kType * NumBin(), NumSys());
	for (int sigma = -3; sigma < 4; sigma += 2)
		sysMatrix[sigma] = Eigen::MatrixXd::Zero(kType * NumBin(),
							 NumSys());

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
			//spline systematic
			if (sysname.find_first_of('_') != sysname.find_last_of('_'))
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
				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number now

				k_err = std::stoi(sysname);
				if (k_err < sysA || k_err >= sysB)
					continue;
				k_err += kScale;

				//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				for (int n = 0; n < NumBin(); ++n)
					sysMatrix[sigma](it * NumBin() + n,
							 k_err - sysA)
						= (hsys->GetBinContent(n+1) - 1);
			}
			else
			{
				TH1D* hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number

				k_err = std::stoi(sysname);
				if (k_err < sysA || k_err >= sysB)
					continue;
				k_err += kScale;

					//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				for (sigma = -3; sigma < 4; sigma += 2)
					for (int n = 0; n < NumBin(); ++n)
						sysMatrix[sigma](it * NumBin() + n,
								k_err - sysA)
							= sigma * (hsys->GetBinContent(n+1) - 1);
			}
		}
		sysf->Close();
	}
}


/*
 * This routine builds the observables as a long vector of NumBin length
 */
Eigen::VectorXd NewChiSquared::ConstructSpectrum(Oscillator *osc)
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


Eigen::VectorXd NewChiSquared::LoadSpectrum(int pt, std::string from)
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

bool NewChiSquared::PointInFile(TFile *f, int pt)
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

int NewChiSquared::NumSys()
{
	return corr.cols();
}

int NewChiSquared::NumBin()
{
	if (_nBin > 0)
		return _nBin;

	std::string reco_file;

	if (!cd->Get("reco_nuE0_nuE0_FHC", reco_file))
		return 0;

	const double *bins;
	Reco *reco = new Reco(reco_file);
	_nBin = reco->BinsY(bins);
	_bins.insert(_bins.end(), bins, bins + _nBin + 1);

	delete reco;

	return _nBin;
}

int NewChiSquared::DOF()
{
	return std::abs(kType * NumBin() - NumSys());
}

//On is the true spectrum, En is the observed spectrum
//return time taken for computation
Eigen::VectorXd NewChiSquared::FitX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	//initialize epsil with zeroes

	Eigen::VectorXd epsil = Eigen::VectorXd::Zero(NumSys());
	Eigen::VectorXd best_eps = epsil;
	double best_x2 = X2(On, En, best_eps);
	
	int tries = 0;
	while (!FindMinimum(On, En, epsil)/*, alpha)*/
		&& tries < maxIteration) {
		++tries;
		double x2 = X2(On, En, epsil);

		if (kVerbosity > 1)
			std::cout << "~~~~~~ No convergence reached ~~~~~~~";

		double dist = (best_eps - epsil).norm();
		double step = x2 - best_x2;
		if (best_x2 > x2) {
			best_x2 = x2;
			best_eps = epsil;
			std::cout << "but new best!";
		}


		if (kVerbosity > 1) {
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
bool NewChiSquared::FindMinimum(const Eigen::VectorXd &On, const Eigen::VectorXd &En, 
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
	//while (delta.norm() > NumSys() * err && c < maxIteration)
	while (std::abs(diff) / DOF() > err && delta.norm() / NumSys() > err) {
		++c;	//counter

		// build hessian and gradient/jacobian together
		// to save computation time
		Eigen::VectorXd jac(NumSys());
		Eigen::MatrixXd hes(NumSys(), NumSys());
		JacobianHessian2(jac, hes, On, En, epsil);


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
		if (x2_ < x2) {	//next x2 is better, update lambda and epsilons
			lambda /= lm_down;
			epsil = nextp;
			x2 = x2_;
		}
		else		//next x2 is worse, nothing changes but lambda
			lambda *= lm_up;

		//if (false)
		if (kVerbosity > 2) {
			std::cout << c << " -> l " << lambda
				  << ",\tstep: " << std::abs(delta.norm() / NumSys())
				  << ", X2: " << x2
				  << " ( " << std::abs(diff / DOF())
				  << " ) " << std::endl;
		}
	}


	if (lambda < 1./lm_0)	//convergence was reached
		return true;
	else			//no convergence, change starting point
		return false;
}


Eigen::VectorXd NewChiSquared::Jacobian(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En, 
				      const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(NumSys());
	Eigen::MatrixXd hes(NumSys(), NumSys());
	JacobianHessian2(jac, hes, On, En, epsil);

	return jac;
}


Eigen::MatrixXd NewChiSquared::Hessian(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En, 
				     const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(NumSys());
	Eigen::MatrixXd hes(NumSys(), NumSys());
	JacobianHessian2(jac, hes, On, En, epsil);

	std::cout << "hessian diagonal\n" << hes.diagonal() << "\n\n";
	return hes+corr;
}

void NewChiSquared::JacobianHessian2(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				    const Eigen::VectorXd &On,
				    const Eigen::VectorXd &En, 
				    const Eigen::VectorXd &epsil)
{
	//corr is inverse of correlation matrix
	jac = corr * epsil;	//gradient/jacobian
	//hes = corr;		//hessian
	hes.setZero();		//hessian

	const Eigen::VectorXd sk  = epsil.head(kScale);

	const Eigen::VectorXd Ep  = Gamma(En, epsil);
	const Eigen::MatrixXd gam = GammaJac(Ep, epsil);

	for (int n = 0; n < kType * NumBin(); ++n) {
		//scale systematic is related to bin number
		int t = n / (kType * NumBin() / kScale);
		double en = Scale(Ep, sk(t), n);

		if (en < 1e-9)		//bad stuff, avoid it
			continue;

		double one_oe = 1 - On(n) / en;
		double en_jac = ScaleJac(Ep, sk(t), n);

		hes(t, t) += On(n) * pow(en_jac / en, 2);
		std::cout << "adding " << n << "\t"
			  << On(n) * pow(en_jac / en, 2)
			  << " += " << hes(t, t) << std::endl;;

		if (std::abs(one_oe) > 1e-9) {
			jac(t)    += one_oe * en_jac;
			hes(t, t) += one_oe * ScaleHes(Ep, sk(t), n);
		}

		//do the rest, icluding mixed term with SK syst
		for (int k = kScale; k < NumSys(); ++k) {
			const Eigen::VectorXd &gk = gam.col(k);
			double gn_jac = Scale(gk, sk(t), n);

			//mixed term with SK error
			hes(k, t) += en_jac * gn_jac * On(n) / en / en;

			//diagonal term
			hes(k, k) += On(n) * pow(gn_jac / en, 2);

			if (std::abs(one_oe) < 1e-9) {
				jac(k)    += one_oe * gn_jac;
				hes(k, t) += one_oe * ScaleJac(gk, sk(t), n);
			}

			for (int j = k + 1; j < NumSys(); ++j) {
				const Eigen::VectorXd &gj = gam.col(j);

				hes(k, j) += gn_jac * Scale(gj, sk(t), n)
					   * On(n) / en / en;

				if (std::abs(one_oe) < 1e-9)
					continue;

				const Eigen::VectorXd gh = GammaJac(gk, epsil, j);
				hes(k, j) += one_oe * Scale(gh, sk(t), n);
			}
		}
	}

	for (int k = 0; k < NumSys(); ++k)
		for (int j = k+1; j < NumSys(); ++j)
			hes(j, k) = hes(k, j);
}



/*
void NewChiSquared::JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				    const Eigen::VectorXd &On,
				    const Eigen::VectorXd &En, 
				    const Eigen::VectorXd &epsil)
{
	//kTpye is also the number of scale errors
	Eigen::VectorXd sk = epsil.head(kScale);

	//modified expected events with systematics
	Eigen::VectorXd Ep = Gamma(En, epsil);
	//Eigen::MatrixXd gk = GammaJac(En, epsil);

	//corr is inverse of correlation matrix
	jac = corr * epsil;	//gradient/jacobian
	hes = corr;		//hessian

	//do SK systematic first (k = 0)
	for (int n = 0; n < kType * NumBin(); ++n) {
		double en = Scale(Ep, sk, n);
		if (en < 1e-9)
			continue;

		int t = n / (kType * NumBin() / kScale);
		double epskn = ScaleJac(Ep, sk, n);

		hes(t, t) += On(n) * pow(epskn / en, 2);

		if (std::abs(en - On(n)) < 1e-9)
			continue;

		jac(t)    += (1 - On(n) / en) * epskn;
		hes(t, t) += (1 - On(n) / en) * ScaleHes(Ep, sk, n);
	}

	//do the rest, icluding mixed term with SK syst
	for (int k = kScale; k < NumSys(); ++k) {
		//gk is jacobian of modified systematic wrt to k
		Eigen::VectorXd gk = GammaJac(En, epsil, k);

		for (int n = 0; n < kType * NumBin(); ++n) {
			double en = Scale(Ep, sk, n);
			if (en < 1e-9)
				continue;

			int t = n / (kType * NumBin() / kScale);
			double gkskn = Scale(gk, sk, n);

			//diagonal term
			hes(k, k) += On(n) * pow(gkskn / en, 2);

			//mixed term with SK error
			hes(k, t) += ScaleJac(Ep, sk, n) * gkskn
				     * On(n) / en / en;

			if (std::abs(en - On(n)) < 1e-9)
				continue;

			//remaining term in mixed term with SK error
			jac(k)    += (1 - On(n) / en) * gkskn;
			hes(k, t) += (1 - On(n) / en) * ScaleJac(gk, sk, n);
		}
		for (int t = 0; t < kType; ++t)
			hes(t, k) = hes(k, t);	//symmetrise


		for (int j = k + 1; j < NumSys(); ++j) {
			Eigen::VectorXd gh = GammaHes(En, epsil, k, j);
			Eigen::VectorXd gj = GammaJac(En, epsil, j);
			for (int n = 0; n < kType * NumBin(); ++n) {
				double en = Scale(Ep, sk, n);
				if (en < 1e-9)
					continue;

				int t = n / (kType * NumBin() / kScale);
				hes(k, j) += Scale(gk, sk, n) * Scale(gj, sk, n)
					   * On(n) / en / en;

				if (std::abs(en - On(n)) < 1e-9)
					continue;
				hes(k, j) += (1 - On(n) / en) * Scale(gh, sk, n);
			}
			hes(j, k) = hes(k, j);
		}
	}
}
*/


// Inverted hessian
Eigen::MatrixXd NewChiSquared::Covariance(const Eigen::VectorXd &On,
					  const Eigen::VectorXd &En,
					  const Eigen::VectorXd &epsil)
{
	//modified expected events with systematics
	Eigen::MatrixXd hes = Hessian(On, En, epsil);

	//return hes.inverse();
	return hes;
}


double NewChiSquared::X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		      const Eigen::VectorXd &epsil)
{
	return ObsX2(On, En, epsil) + SysX2(epsil);
}


// epsil is an array with all sigmas including SK energy scale
// SK energy scale is the first one
double NewChiSquared::ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			    const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	Eigen::VectorXd gam = Gamma(En, epsil);
	Eigen::VectorXd sk = epsil.head(kScale);


	std::vector<double> chi2n = ObsX2n(On, En, epsil);
	double chi2 = 0;
	for (int n = 0; n < chi2n.size(); ++n)
		chi2 += chi2n[n];

	return chi2;
}


std::vector<double> NewChiSquared::ObsX2n(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	Eigen::VectorXd gam = Gamma(En, epsil);
	Eigen::VectorXd sk = epsil.head(kScale);

	std::vector<double> chi2(kType * NumBin());
	for (int n = 0; n < kType * NumBin(); ++n) {
		int t = n / (kType * NumBin() / kScale);
		double en = Scale(gam, sk(t), n);
		if (On(n) > 1e-15) {
			en = std::max(1e-8, en);
			chi2[n] = 2 * On(n) * ( en/On(n) + log(On(n)/en) - 1);
		}
		else
			chi2[n] = 2 * (en - On(n));
	}

	return chi2;
}



double NewChiSquared::SysX2(const Eigen::VectorXd &epsil)
{
	return epsil.transpose() * corr * epsil;
}


// this return the spectrum modified with the systematics
Eigen::VectorXd NewChiSquared::Gamma(const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd gam = En;
	for (int k = kScale; k < NumSys(); ++k)
		for (int n = 0; n < kType * NumBin(); ++n)
			gam(n) *= 1 + F(k, n, epsil(k));

	return gam;
}


// this is the derivative of modified spectrum by gamma at given k
Eigen::MatrixXd NewChiSquared::GammaJac(const Eigen::VectorXd &En,
				        const Eigen::VectorXd &epsil)
{
	Eigen::MatrixXd gam(En.size(), NumSys());
	gam.leftCols(kScale).setZero();
	for (int k = kScale; k < NumSys(); ++k)
		gam.col(k) = GammaJac(En, epsil, k);

	return gam;
}


// derive by k-th systematic
// En is already derived once, so this can be done recursively
Eigen::VectorXd NewChiSquared::GammaJac(const Eigen::VectorXd &En,
				        const Eigen::VectorXd &epsil, int k)
{
	Eigen::VectorXd gam = En;
	for (int n = 0; n < gam.size(); ++n)
		gam(n) *= Fp(k, n, epsil(k)) / (1 + F(k, n, epsil(k)));

	return gam;
}


//F has the epsilon in it
double NewChiSquared::F(int k, int n, double eij)
{
	double dl, du;

	if (eij < -1) {
		dl = -3; du = -1;
	}
	else if (eij >= -1 && eij < 0) {
		dl = -1; du = 0;
	}
	else if (eij >= 0 && eij < 1) {
		dl = 0; du = 1;
	}
	else if (eij >= 1) {
		dl = 1; du = 3;
	}

	return Fp(k, n, dl, du) * (eij - dl) + sysMatrix[dl](n, k);
}

//Fp does not have the epsilon in it
double NewChiSquared::Fp(int k, int n, double eij)
{
	double dl, du;

	if (eij < -1) {
		dl = -3; du = -1;
	}
	else if (eij >= -1 && eij < 0) {
		dl = -1; du = 0;
	}
	else if (eij >= 0 && eij < 1) {
		dl = 0; du = 1;
	}
	else if (eij >= 1) {
		dl = 1; du = 3;
	}

	return Fp(k, n, dl, du);
}

double NewChiSquared::Fp(int k, int n, double dl, double du)
{
	return (sysMatrix[du](n, k) - sysMatrix[dl](n, k)) / (du - dl);
}



// binary search for starting bin
// makes the computation slightly faster
int NewChiSquared::StartingBin(int n, double scale)
{
	double b0_n = _bins[n];
	double b1_n = _bins[n + 1];

	int m0 = 0;
	int m1 = NumBin();
	int mm = (m0 + m1) / 2;
	while (scale * _bins[mm] > b0_n || scale * _bins[mm+1] < b0_n) {
		if (b0_n <= scale * _bins[mm + 1]) {
			m0 = m0;
			m1 = mm;
		}
		else if (b0_n >= scale * _bins[mm]) {
			m0 = mm;
			m1 = m1;
		}
		else
			break;
		mm = (m0 + m1) / 2;
	}
	return mm;
}

// Apply energy scale dilation
// The dilation does not depend on En, therefore 
// En can be anything
//double NewChiSquared::Scale(FactorFn factor, const Eigen::VectorXd &En,
//			const Eigen::VectorXd &sigma, int t)
double NewChiSquared::Scale(FactorFn factor, const Eigen::VectorXd &En,
			double sigma, int t)
{
	if (t < 0 || t > kType * NumBin())
		return 0;

	int n = t % NumBin();
	t /= NumBin();

	if (t >= kType)
		return 0;

	double SKerror = 0;
	for (is = scale.begin(); is != scale.end(); ++is) 
		if (type[t].find(is->first) != std::string::npos) {
			SKerror = is->second;
			break;
		}
	double scale = 1 + sigma * SKerror;

	// unscaled/original bins
	double b0_n = _bins[n];
	double b1_n = _bins[n + 1];

	//binary search for starting bin
	int mm = StartingBin(n, scale);

	double ret = 0;
	// Loop over unscaled/original edges
	for (int m = mm; m < NumBin(); ++m) {

		// scaled bins
		double b0_m = _bins[m];
		double b1_m = _bins[m + 1];

		//stop computation because there is no overlap
		if (scale * b0_m > b1_n)
			break;

		double f = (b1_n - b0_n + scale * (b1_m - b0_m)
			  - std::abs(b0_n - scale * b0_m)
			  - std::abs(b1_n - scale * b1_m)) / 2.;

		if (f < 0)
			continue;

		int s0 = b0_n - scale * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - scale * b1_m < 0 ? -1 : 1;
		double fd = (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;
		double ss = f > 0 ? 1 : 0.5;	//continuity factor

		std::vector<double> terms = {SKerror, scale, f, fd, b1_m - b0_m};
		ret += En(m + t * NumBin()) * ss * (this->*factor)(terms);

		//normal
		//ret += En(m + t * NumBin()) * f / scale / (b1_m - b0_m) / 2.;
		//factor
		//return f / scale / db / 2.;

		//jacobian
		//ret += En(m + t * NumBin()) * SKerror / scale / 2.
		//	    / (b1_m - b0_m) * (fd - f / scale);
		//factor
		//return SKerror / scale / 2. / db * (fd - f / scale)

		//hessian
		//ret += En(t) * pow(SKerror / scale, 2)
		//	/ (b1_m - b0_m) * (f / scale - fd);
		//factor
		//return pow(SKerror / scale, 2) / db * (f / scale - fd)
	}

	if (factor == &NewChiSquared::ScaleNor
	 && n == NumBin() - 1 &&  scale > 1) { //to keep same integral
		double b0_n = _bins[n];
		double b1_n = _bins[n + 1];

		// "1-" because this is the bit that is above
		// the maximum (e.g. 30 GeV for 1Rmu)
		// We've already counted the other bit
		double f = 1 - (b1_n - scale * b0_n)
			   / scale / (b1_n - b0_n);
		ret += En(n + t * NumBin()) * f;
	}

	return ret;
}


double NewChiSquared::Scale(const Eigen::VectorXd &En, double sigma, int t)
{
	return Scale(&NewChiSquared::ScaleNor, En, sigma, t);
}


double NewChiSquared::ScaleJac(const Eigen::VectorXd &En, double sigma, int t)
{
	return Scale(&NewChiSquared::ScaleJac, En, sigma, t);
}


double NewChiSquared::ScaleHes(const Eigen::VectorXd &En, double sigma, int t)
{
	return Scale(&NewChiSquared::ScaleHes, En, sigma, t);
}


double NewChiSquared::ScaleNor(const std::vector<double> &term)
{
	//return f / scale / db;
	return term[2] / term[1] / term[4];
}

double NewChiSquared::ScaleJac(const std::vector<double> &term)
{
	//return SKerror / scale / db * (fd - f / scale)
	return term[0] / term[1] / term[4] * (term[3] - term[2] / term[1]);
}

double NewChiSquared::ScaleHes(const std::vector<double> &term)
{
	//return pow(SKerror / scale, 2) / db * (f / scale - fd)
	return 2 * pow(term[0] / term[1], 2) / term[4]
		 * (term[2] / term[1] - term[3]);
}



/*
// returns the derivative with respect to sigma
// can be used with non diagonal hessian terms involving sigma
// just replacing En with the right vector
double NewChiSquared::ScaleJac(const Eigen::VectorXd &En,
			       const Eigen::VectorXd &sigma, int t)
{
	if (t < 0 || t > kType * NumBin())
		return 0;

	int n = t % NumBin();
	t /= (kType * NumBin() / kScale);

	if (t >= kType)
		return 0;

	double SKerror = 0;
	std::map<std::string, double>::iterator is;
	for (is = scale.begin(); is != scale.end(); ++is) 
		if (type[t].find(is->first) != std::string::npos) {
			SKerror = is->second;
			break;
		}
	double scale = 1 + sigma(t) * SKerror;

	// unscaled/original bins
	double b0_n = _bins[n];
	double b1_n = _bins[n + 1];

	double ret = 0;
	// Loop over unscaled/original edges
	for (int m = 0; m < NumBin(); ++m) {

		// scaled bins
		double b0_m = _bins[m];
		double b1_m = _bins[m + 1];

		double f = b1_n - b0_n + scale * (b1_m - b0_m)
			   - std::abs(b0_n - scale * b0_m)
			   - std::abs(b1_n - scale * b1_m);

		if (f <= 0)
			continue;

		int s0 = b0_n - scale * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - scale * b1_m < 0 ? -1 : 1;
		double fd = b1_m - b0_m + s0 * b0_m + s1 * b1_m;

		ret += En(m + t * NumBin()) * SKerror / scale / 2.
			    / (b1_m - b0_m) * (fd - f / scale);
		if (!std::isfinite(ret))
			std::cout << "ff " << En(m + t * NumBin()) << "\t" << (b1_m - b0_m) << "\t" << fd << "\t" <<  f << std::endl;
	}

	return ret;
}

// return double derivative on SK error
double NewChiSquared::ScaleHes(const Eigen::VectorXd &En,
			       const Eigen::VectorXd &sigma, int t)
{
	if (t < 0 || t > kType * NumBin())
		return 0;

	int n = t % NumBin();
	t /= (kType * NumBin() / kScale);

	if (t >= kType)
		return 0;

	double SKerror = 0;
	std::map<std::string, double>::iterator is;
	for (is = scale.begin(); is != scale.end(); ++is) 
		if (type[t].find(is->first) != std::string::npos) {
			SKerror = is->second;
			break;
		}
	double scale = 1 + sigma(t) * SKerror;

	// unscaled/original bins
	double b0_n = _bins[n];
	double b1_n = _bins[n + 1];

	double ret = 0;
	// Loop over unscaled/original edges
	for (int m = 0; m < NumBin(); ++m) {

		// scaled bins
		double b0_m = _bins[m];
		double b1_m = _bins[m + 1];

		double f = b1_n - b0_n + scale * (b1_m - b0_m)
			- std::abs(b0_n - scale * b0_m)
			- std::abs(b1_n - scale * b1_m);

		if (f <= 0)
			continue;

		int s0 = b0_n - scale * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - scale * b1_m < 0 ? -1 : 1;
		double fd = b1_m - b0_m + s0 * b0_m + s1 * b1_m;

		ret += En(t) * pow(SKerror / scale, 2)
		       / (b1_m - b0_m) * (f / scale - fd);
	}

	return ret;
}
*/
