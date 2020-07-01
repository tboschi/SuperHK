/* SK energy scale smart handling
 *
 * It was decided that bins with zero events should be not inclded in chi^2 calculation
 * especially due to the energy scale problem
 * If events migrate to previously-empty bins, the chi^2 becomes ill-defenied.
 *
 * For this reason, it should be safer (and faster too!!) to only store and handle nonzero bins
 * 1Re-like hists and 1Rmu-like hists have therefore a different binning range
 */

#include "ChiSquared.h"

ChiSquared::ChiSquared(CardDealer *card, int nbins) :
	cd(card),
	_nBin(nbins),
	_nSys(-1)
{
	_mode = std::vector<std::string>
		({"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0", "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"});
	_chan = std::vector<std::string>
		({"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"});
	_horn = std::vector<std::string>
		({"FHC", "RHC"});
	_oscf = std::vector<std::string>
		({"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"});
	 _fin = std::vector<Nu> ({Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb});
	_fout = std::vector<Nu> ({Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb});

	if (!cd->Get("sample", _type)) {
		if (nbins < 0)
			_type = std::vector<std::string> ({"E_FHC", "E_RHC", "M_FHC", "M_RHC"});
		else
			_type = std::vector<std::string>(1);
	}
	else	//sample is defined

	// maybe one scale error?	
	if (!cd->Get("scale_", _scale)) {
		double err;
		_scale.clear();
		if (!cd->Get("scale", err))
			err = 0.024;

		for (const std::string &it : _type) {
			_scale[it] = err;
			_type_scale[it] = 0;
		}
	}
	else {	// two scale_E and scale_M
		// scale contains "E" and "M"	
		std::map<std::string, double> err;
		for (const std::string &it : _type) {
			if (it.find_first_of('E') != std::string::npos) {
				if (_scale.count("E"))
					err[it] = _scale["E"];
				else {
					if (kVerbosity)
						std::cout << "No scale_E in card"
							  << cd->CardName() << std::endl;
					err[it] = 0.024;
				}
				_type_scale[it] = 0;
			}
			else if (it.find_first_of('M') != std::string::npos) {
				if (_scale.count("M"))
					err[it] = _scale["M"];
				else {
					if (kVerbosity)
						std::cout << "No scale_M in card"
							  << cd->CardName() << std::endl;
					err[it] = 0.024;
				}
				_type_scale[it] = 1;
			}
		}
		// reassign new map
		_scale = err;
	}
	_nScale = _scale.size();	//of size 1, 2, or 4


	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "types: ";
		for (int t = 0; t < _type.size(); ++t)
			std::cout << "\t" << _type[t];
		std::cout << std::endl;
		std::cout << "scales: " << _nScale << std::endl;
		for (const auto &it : _scale)
			std::cout << "\t" << it.first << " -> " << it.second << std::endl;
	}

	maxIteration = 10;
	err = 1e-9;

	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	spectrumFile = 0;
	input = "input";
}


int ChiSquared::NumSys()
{
	if (_nSys > 0)
		return _nSys;

	_nSys = corr.cols();
	return _nSys;
}


int ChiSquared::NumBin() {
	// compute the number of nonempty bins
	// if not defined, the routine will do this check
	// and creates lower and upper bin pairs in _limits
	//
	if (_nBin > 0)	// all number of bins (compressed)
		return _nBin;

	DefineBinning();
	return _nBin;
}

void ChiSquared::DefineBinning() {
	// extract bining from spectra
	// binning is compressed, such that nonempty bins are not stored

	// unoscillated spectrum, just for counting nonempty bins
	std::map<std::string, TH1D*> spectra = BuildSpectrum();

	_nBin = 0;
	_allBin = 0;
	_limits.clear();
	_bins.clear();
	_global.clear();

	for (const auto &it : spectra) {
		int b0 = std::max(it.second->FindFirstBinAbove() - 1, 0);
		int b1 = std::max(it.second->FindLastBinAbove(), it.second->GetNbinsX());

		// store first and last bin
		// refers to TH1D
		_limits[it.first] = std::pair<int, int>(b0, b1);

		// store binning as vector in map
		const double *bins = it.second->GetXaxis()
					->GetXbins()->GetArray();
		_bins[it.first].assign(bins, bins + it.second->GetNbinsX()
						  + 1);

		// refers to Eigen
		std::vector<int> binmap(b1-b0);
		std::iota(binmap.begin(), binmap.end(), b0);
		_global.insert(_global.end(), binmap.begin(), binmap.end());

		// count number of total bins
		_nBin += b1 - b0;
		_allBin += it.second->GetNbinsX();

		delete it.second;
	}

	if (kVerbosity)
		std::cout << "Number of nonzero bins: "
			  << _allBin << " / " << _nBin << std::endl;
}

int ChiSquared::DOF() {
	return std::abs(_nBin - _nSys);
}


void ChiSquared::Init()
{
	DefineBinning();

	LoadCorrelation();
	LoadSystematics();
}

void ChiSquared::LoadCorrelation()
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
	//last _nScale entries are for energy scale
	corr = Eigen::MatrixXd::Zero(_nScale + sysB - sysA, _nScale + sysB - sysA);
	corr.bottomRightCorner(_nScale, _nScale) = Eigen::MatrixXd::Identity(_nScale, _nScale);

	for (int r = sysA; r < sysB; ++r)
		for (int c = sysA; c < sysB; ++c)
			corr(r + - sysA, c + - sysA) = (*cmat)(r, c);
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
void ChiSquared::LoadSystematics()
{
	sysMatrix.clear();
	sysMatrix[0] = Eigen::MatrixXd::Zero(_nBin, _nSys);
	for (int sigma = -3; sigma < 4; sigma += 2)
		sysMatrix[sigma] = Eigen::MatrixXd::Zero(_nBin, _nSys);

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	if (!cd->Get("syst_last", sysB))
		sysB = _nSys;

	int off = 0;	// offset for global bin
	for (const std::string &it : _type) {
		// it is sample type name
		std::string fileName;
		if (!cd->Get("systematic_" + it, fileName))
			continue;

		if (kVerbosity > 1)
			std::cout << "opening file " << fileName << std::endl;
		TFile *sysf = new TFile(fileName.c_str());
		TIter next(sysf->GetListOfKeys());
		TKey *k;

		// bin limits
		auto lims = _limits[it];
		while (k = static_cast<TKey*>(next()))
		{
			std::string sysname = k->GetName();

			//if (skip_syst.find(sysname) != skip_sys.end()) {
			//	if (kVerbosity)
			//		std::cout << "Not accepting " << sysname << std::endl;
			//	continue;	// to be implemented
			//}
			//elif (kVerbosity)
			//	std::cout << "Accepting " << sysname << std::endl;

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
				//k_err += _nScale;

				//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				for (int n = lims.first, j = 0; n < lims.second; ++n, ++j)
					sysMatrix[sigma](off + j, k_err - sysA)
						= hsys->GetBinContent(n+1) - 1;
			}
			else
			{
				TH1D* hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number

				k_err = std::stoi(sysname);
				if (k_err < sysA || k_err >= sysB)
					continue;
				//k_err += _nScale;

					//for (int n = 0; n < hsys->GetNbinsX(); ++n)
				//for (int n = 0; n < NumBin(); ++n)
				for (int n = lims.first, j = 0; n < lims.second; ++n, ++j)
					for (sigma = -3; sigma < 4; sigma += 2)
						sysMatrix[sigma](off + j, k_err - sysA)
							= sigma * (hsys->GetBinContent(n+1) - 1);
			}
		}
		off += lims.second - lims.first;
		sysf->Close();
	}
}


Eigen::VectorXd ChiSquared::ConstructSpectrum(Oscillator *osc) {
	// BuildSpectrum and then collates everything on a Eigen::Vector

	std::map<std::string, TH1D*> spectra = BuildSpectrum(osc);

	double stats;	//for scaling
	if (!cd->Get("stats", stats))
		stats = 1.0;

	Eigen::VectorXd vect(_nBin);

	int j = 0;	//add offset
	for (const auto &is : spectra) {
		auto lims = _limits[is.first];
		for (int n = lims.first; n < lims.second; ++n, ++j) 
			vect(j) = is.second->GetBinContent(n+1) * stats;

		delete is.second;
	}

	return vect;
}

std::map<std::string, TH1D*> ChiSquared::BuildSpectrum(Oscillator *osc) {
	// build the observables as TH1D objects
	// and return map with histograms for each type
	// if no oscillator is passed (default = 0), spectra are not oscillated
	//
	std::string reco_file;

	TF1* f1 = new TF1("const", "[0]", 0, 1);
	f1->SetParameter(0, 1.0);

	std::map<std::string, TH1D*> spectra;

	for (const std::string &ih : _horn) {	//FHC, RHC
		int cm = 0;	// mode counter
		for (const std::string &im : _mode) {	//nuE->nuE, nuM->nuE, nuM->nuM
			// Reco object combines reconstruction matrices to for prediction
			cd->Get("reco_" + im + "_" + ih, reco_file);
			Reco *reco = new Reco(reco_file);

			const double *bins;
			int nBin = reco->BinsX(bins);
			f1->SetRange(bins[0], bins[nBin-1]);

			// flux contains the weight of spectrum
			std::string fname = ih + "_" + im;
			TH1D* flux = new TH1D(fname.c_str(), fname.c_str(), nBin, bins);
			flux->SetDirectory(0);
			flux->Add(f1);	//filled with ones

			if (osc)	
				osc->Oscillate(_fin[cm], _fout[cm], flux);

			for (const std::string &ic : _chan) {	//CCQE, CCnQE, NC x E, M
				if (ic.find("NC") == std::string::npos)	//it is not a NC
					reco->Scale(ic, flux);
				else if (im == "nuM0_nuE0" || im == "nuMB_nuEB") //and no events for
					reco->Scale(ic, 0.0);			     //these two
				//else no need to scale for NC event
				//reco->Scale(chan[ic], flux);

				TH1D *py = reco->Project(ic, 'y');
				if (py) {
					std::string hname = ic.substr(0, ic.find_first_of('_')) + "_" + ih;

					if (spectra.count(hname) && spectra[hname])
						spectra[hname]->Add(py);
					else
						spectra[hname] = static_cast<TH1D*>(py->Clone());
				}
			}
			delete flux;
			delete reco;
		}
		++cm;
	}

	delete f1;

	return spectra;
}

// Atmospheric function
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


	Eigen::VectorXd vect(_nBin);

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



//On is the true spectrum, En is the observed spectrum
//return time taken for computation
Eigen::VectorXd ChiSquared::FitX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	//initialize epsil with zeroes

	Eigen::VectorXd epsil = Eigen::VectorXd::Zero(_nSys);
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

	Eigen::VectorXd delta = Eigen::VectorXd::Constant(_nSys, 1);

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


Eigen::VectorXd ChiSquared::Jacobian(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En, 
				      const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(NumSys());
	Eigen::MatrixXd hes(NumSys(), NumSys());
	JacobianHessian2(jac, hes, On, En, epsil);

	return jac;
}


Eigen::MatrixXd ChiSquared::Hessian(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En, 
				     const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(NumSys());
	Eigen::MatrixXd hes(NumSys(), NumSys());
	JacobianHessian2(jac, hes, On, En, epsil);

	std::cout << "hessian diagonal\n" << hes.diagonal() << "\n\n";
	return hes+corr;
}

void ChiSquared::JacobianHessian2(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				    const Eigen::VectorXd &On,
				    const Eigen::VectorXd &En, 
				    const Eigen::VectorXd &epsil)
{
	//corr is inverse of correlation matrix
	jac = corr * epsil;	//gradient/jacobian
	//hes = corr;		//hessian
	hes.setZero();		//hessian

	const Eigen::VectorXd sk  = epsil.tail(_nScale);

	const Eigen::VectorXd Ep  = Gamma(En, epsil);
	const Eigen::MatrixXd gam = GammaJac(Ep, epsil);

	for (int n = 0; n < _nBin; ++n) {
		// scale systematic is related to bin number
		int t = _nSys - _nScale + _type_scale[TypeFromBin(n)];

		double en = Scale(Ep, sk, n);	// en can be empty!
		if (en < 1e-9)			//bad stuff, avoid it
			continue;

		double one_oe = 1 - On(n) / en;
		double en_jac = ScaleJac(Ep, sk, n);

		// hessian of scale error first
		hes(t, t) += On(n) * pow(en_jac / en, 2);
		//std::cout << "adding " << n << "\t"
			  //<< On(n) * pow(en_jac / en, 2)
			  //<< " += " << hes(t, t) << std::endl;;

		if (std::abs(one_oe) > 1e-9) {
			jac(t)    += one_oe * en_jac;
			hes(t, t) += one_oe * ScaleHes(Ep, sk, n);
		}

		//do the rest, including mixed term with SK syst
		for (int k = 0; k < _nSys - _nScale; ++k) {
			const Eigen::VectorXd &gk = gam.col(k);
			double gn_jac = Scale(gk, sk, n);

			//mixed term with SK error
			hes(k, t) += en_jac * gn_jac * On(n) / en / en;

			//diagonal term
			hes(k, k) += On(n) * pow(gn_jac / en, 2);

			if (std::abs(one_oe) < 1e-9) {
				jac(k)    += one_oe * gn_jac;
				hes(k, t) += one_oe * ScaleJac(gk, sk, n);
			}

			for (int j = k + 1; j < _nSys - _nScale; ++j) {
				const Eigen::VectorXd &gj = gam.col(j);

				hes(k, j) += gn_jac * Scale(gj, sk, n)
					   * On(n) / en / en;

				if (std::abs(one_oe) < 1e-9)
					continue;

				const Eigen::VectorXd gh = GammaJac(gk, epsil, j);
				hes(k, j) += one_oe * Scale(gh, sk, n);
			}
		}
	}

	for (int k = 0; k < _nSys - _nScale; ++k)
		for (int j = k+1; j < _nSys - _nScale; ++j)
			hes(j, k) = hes(k, j);
}



/*
void ChiSquared::JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				    const Eigen::VectorXd &On,
				    const Eigen::VectorXd &En, 
				    const Eigen::VectorXd &epsil)
{
	//kTpye is also the number of scale errors
	Eigen::VectorXd sk = epsil.head(_nScale);

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

		int t = n / (kType * NumBin() / _nScale);
		double epskn = ScaleJac(Ep, sk, n);

		hes(t, t) += On(n) * pow(epskn / en, 2);

		if (std::abs(en - On(n)) < 1e-9)
			continue;

		jac(t)    += (1 - On(n) / en) * epskn;
		hes(t, t) += (1 - On(n) / en) * ScaleHes(Ep, sk, n);
	}

	//do the rest, icluding mixed term with SK syst
	for (int k = _nScale; k < NumSys(); ++k) {
		//gk is jacobian of modified systematic wrt to k
		Eigen::VectorXd gk = GammaJac(En, epsil, k);

		for (int n = 0; n < kType * NumBin(); ++n) {
			double en = Scale(Ep, sk, n);
			if (en < 1e-9)
				continue;

			int t = n / (kType * NumBin() / _nScale);
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

				int t = n / (kType * NumBin() / _nScale);
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
Eigen::MatrixXd ChiSquared::Covariance(const Eigen::VectorXd &On,
					  const Eigen::VectorXd &En,
					  const Eigen::VectorXd &epsil)
{
	//modified expected events with systematics
	Eigen::MatrixXd hes = Hessian(On, En, epsil);

	//return hes.inverse();
	return hes;
}


double ChiSquared::X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		      const Eigen::VectorXd &epsil)
{
	return ObsX2(On, En, epsil) + SysX2(epsil);
}


// epsil is an array with all sigmas including SK energy scale
// sum up all bin contributions to the chi2
double ChiSquared::ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			    const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	std::vector<double> chi2n = ObsX2n(On, En, epsil);
	return std::accumulate(chi2n.begin(), chi2n.end(), 0);
}


// return statistics X2 as a vector
// so X2 contribution to each bin
std::vector<double> ChiSquared::ObsX2n(const Eigen::VectorXd &On,
				       const Eigen::VectorXd &En,
				       const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	Eigen::VectorXd gam = Gamma(En, epsil);
	Eigen::VectorXd sk = epsil.tail(_nScale);

	std::vector<double> chi2(En.size());
	for (int n = 0; n < En.size(); ++n) {
		double en = std::max(Scale(gam, sk, n), 1e-8);
		chi2[n] = 2 * On(n) * (en / On(n) - log(en / On(n)) - 1);
	}

	return chi2;
}


// return systematic X2
double ChiSquared::SysX2(const Eigen::VectorXd &epsil) {
	return epsil.transpose() * corr * epsil;
}


// this return the spectrum modified with the systematics
Eigen::VectorXd ChiSquared::Gamma(const Eigen::VectorXd &En,
				  const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd gam = En;
	for (int k = 0; k < _nSys - _nScale; ++k)
		for (int n = 0; n < En.size(); ++n)
			gam(n) *= 1 + F(k, n, epsil(k));

	return gam;
}


// this is the derivative of modified spectrum by gamma at given k
Eigen::MatrixXd ChiSquared::GammaJac(const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil)
{
	Eigen::MatrixXd gam(_nBin, _nSys);	// should remove _nScale?
	gam.rightCols(_nScale).setZero();
	for (int k = 0; k < _nSys - _nScale; ++k)
		gam.col(k) = GammaJac(En, epsil, k);

	return gam;
}


// derive by k-th systematic
// En is already derived once, so this can be done recursively
Eigen::VectorXd ChiSquared::GammaJac(const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil, int k)
{
	Eigen::VectorXd gam = En;
	for (int n = 0; n < En.size(); ++n)
		gam(n) *= Fp(k, n, epsil(k)) / (1 + F(k, n, epsil(k)));

	return gam;
}


//F has the epsilon in it
double ChiSquared::F(int k, int n, double eij)
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
double ChiSquared::Fp(int k, int n, double eij)
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

double ChiSquared::Fp(int k, int n, double dl, double du)
{
	return (sysMatrix[du](n, k) - sysMatrix[dl](n, k)) / (du - dl);
}



// binary search for starting bin
// makes the computation slightly faster
/*
int ChiSquared::StartingBin(int n, double scale)
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
*/

// Apply energy scale dilation
// FactorFn is a function which changes if calculating energy scale, jacobian or hessian
// The dilation does not depend on En, therefore 
// En can be anything
//double ChiSquared::Scale(FactorFn factor, const Eigen::VectorXd &En,
//			const Eigen::VectorXd &sigma, int t)
/*
double ChiSquared::Scale2(FactorFn factor, const Eigen::VectorXd &En,
			double sigma, int t)
{
	if (t < 0 || t > NumBin())
		return 0;

	int n = t % NumBin();
	t /= NumBin();

	if (t >= kType)
		return 0;

	double SKerror = 0;
	for (is = _scale.begin(); is != _scale.end(); ++is) 
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

	if (factor == &ChiSquared::ScaleNor
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
*/

std::string ChiSquared::TypeFromBin(int n) {
	// check if is the same type
	if (n >= _limits[_tlim].first && n < _limits[_tlim].second)
		return _tlim;

	// return sample type given bin
	for (const auto &it : _limits)
		if (it.first != _tlim
		    && n >= it.second.first
		    && n < it.second.second) {
			_tlim = it.first;
			return _tlim;
		}

	return "";
}

double ChiSquared::Scale(FactorFn factor,
			 const Eigen::VectorXd &En,
			 const Eigen::VectorXd &sk,
			 int n)	// relative bin
{
	std::string it = TypeFromBin(n);

	double scale_err = _scale[it];	// this is Error value
	double shift = 1 + sk(_type_scale[it]) * scale_err;

	auto lims = _limits[it];
	int off = _global[n] - n;
	n = _global[n];	// global bin

	// unscaled/original bins
	double b0_n = _bins[it][n];
	double b1_n = _bins[it][n + 1];

	//binary search for starting bin
	//int mm = StartingBin(n, shift);

	double ret = 0;
	// Loop over unscaled/original edges only nonempty bins
	for (int m = lims.first; m < lims.second; ++m) {

		// scaled bins
		double b0_m = _bins[it][m];
		double b1_m = _bins[it][m + 1];

		//continue/break computation because there is no overlap
		if (shift * b1_m < b0_n)
			continue;
		else if (shift * b0_m > b1_n)
			break;

		double f = (b1_n - b0_n + shift * (b1_m - b0_m)
			  - std::abs(b0_n - shift * b0_m)
			  - std::abs(b1_n - shift * b1_m)) / 2.;

		if (f < 0)
			continue;

		int s0 = b0_n - shift * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - shift * b1_m < 0 ? -1 : 1;
		double fd = (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;
		double ss = f > 0 ? 1 : 0.5;	//continuity factor

		std::vector<double> terms = {scale_err, shift, f, fd, b1_m - b0_m};
		ret += En(m - off) * ss * (this->*factor)(terms);
	}

	return ret;
}


// call Scale Energy function
double ChiSquared::Scale(const Eigen::VectorXd &En,
			 const Eigen::VectorXd &sk,
			 int n)
{
	return Scale(&ChiSquared::ScaleNor, En, sk, n);
}

// call jacobian of Scale Energy function
double ChiSquared::ScaleJac(const Eigen::VectorXd &En,
			    const Eigen::VectorXd &sk,
			    int n)
{
	return Scale(&ChiSquared::ScaleJac, En, sk, n);
}


// call hessian of Scale Energy function
double ChiSquared::ScaleHes(const Eigen::VectorXd &En,
			    const Eigen::VectorXd &sk,
			    int n)
{
	return Scale(&ChiSquared::ScaleHes, En, sk, n);
}


double ChiSquared::ScaleNor(const std::vector<double> &term)
{
	//return f / shift / db;
	return term[2] / term[1] / term[4];
}

double ChiSquared::ScaleJac(const std::vector<double> &term)
{
	//return scale_err / shift / db * (fd - f / shift)
	return term[0] / term[1] / term[4] * (term[3] - term[2] / term[1]);
}

double ChiSquared::ScaleHes(const std::vector<double> &term)
{
	//return pow(scale_err / shift, 2) / db * (f / shift - fd)
	return 2 * pow(term[0] / term[1], 2) / term[4]
		 * (term[2] / term[1] - term[3]);
}


/*
// returns the derivative with respect to sigma
// can be used with non diagonal hessian terms involving sigma
// just replacing En with the right vector
double ChiSquared::ScaleJac(const Eigen::VectorXd &En,
			       const Eigen::VectorXd &sigma, int t)
{
	if (t < 0 || t > kType * NumBin())
		return 0;

	int n = t % NumBin();
	t /= (kType * NumBin() / _nScale);

	if (t >= kType)
		return 0;

	double scale_err = 0;
	std::map<std::string, double>::iterator is;
	for (is = scale.begin(); is != scale.end(); ++is) 
		if (type[t].find(is->first) != std::string::npos) {
			scale_err = is->second;
			break;
		}
	double scale = 1 + sigma(t) * scale_err;

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

		ret += En(m + t * NumBin()) * scale_err / scale / 2.
			    / (b1_m - b0_m) * (fd - f / scale);
		if (!std::isfinite(ret))
			std::cout << "ff " << En(m + t * NumBin()) << "\t" << (b1_m - b0_m) << "\t" << fd << "\t" <<  f << std::endl;
	}

	return ret;
}

// return double derivative on SK error
double ChiSquared::ScaleHes(const Eigen::VectorXd &En,
			       const Eigen::VectorXd &sigma, int t)
{
	if (t < 0 || t > kType * NumBin())
		return 0;

	int n = t % NumBin();
	t /= (kType * NumBin() / _nScale);

	if (t >= kType)
		return 0;

	double scale_err = 0;
	std::map<std::string, double>::iterator is;
	for (is = scale.begin(); is != scale.end(); ++is) 
		if (type[t].find(is->first) != std::string::npos) {
			scale_err = is->second;
			break;
		}
	double scale = 1 + sigma(t) * scale_err;

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

		ret += En(t) * pow(scale_err / scale, 2)
		       / (b1_m - b0_m) * (f / scale - fd);
	}

	return ret;
}
*/
