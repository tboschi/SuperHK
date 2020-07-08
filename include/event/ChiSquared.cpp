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
	_mode = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0",
		 "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	_chan = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	_horn = {"FHC", "RHC"};
	_oscf = {"oscEE_0", "oscMM_0", "oscME_0", "oscEE_B", "oscMM_B", "oscME_B"};
	 _fin = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	_fout = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};

	if (!cd->Get("sample", _type)) {
		if (nbins < 0)
			_type = {"E_FHC", "E_RHC", "M_FHC", "M_RHC"};
		else
			_type = {""};
	}
	else	//sample is defined

	// maybe one scale error?	
	if (!cd->Get("scale_", _scale)) {
		double err;
		_scale.clear();
		if (cd->Get("scale", err)) {
			_nScale = 1;
			for (const std::string &it : _type) {
				_scale[it] = err;
				_type_scale[it] = 0;
			}
		}
		else
			_nScale = 0;
	}
	else {	// two scale_E and scale_M
		// scale contains "E" and "M"	
		std::map<std::string, double> err;
		for (const std::string &it : _type) {
			//E_FHC or E_RHC
			if (it.find_first_of('E') != std::string::npos) {
				if (_scale.count("E"))
					err[it] = _scale["E"];
				else {
					if (kVerbosity)
						std::cout << "No scale_E in card"
							  << cd->CardName() << std::endl;
					//default to
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
					// default to
					err[it] = 0.024;
				}
				_type_scale[it] = 1;
			}
		}
		_nScale = _scale.size();	//of size 1, 2, or 4
		// reassign new map
		_scale = err;
	}


	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "types: ";
		for (int t = 0; t < _type.size(); ++t)
			std::cout << "\t" << _type[t];
		std::cout << std::endl;
		std::cout << "number of scale systematics " << _nScale << std::endl;
		for (const auto &it : _scale)
			std::cout << "\t" << it.first << " -> " << it.second << std::endl;
	}

	maxIteration = 10;
	fitErr = 1e-6;

	if (!cd->Get("lm_0", lm_0))
		lm_0 = 1;		//default value
	if (!cd->Get("lm_up", lm_up))
		lm_up = 5;		//default value
	if (!cd->Get("lm_down", lm_down))
		lm_down = 10;		//default value


	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;

	spectrumFile = 0;
	input = "input";

	Init();
}


int ChiSquared::NumSys()
{
	if (_nSys < 0)
		std::cerr << "WARNING - ChiSquared: systematics not set!" << std::endl;

	return _nSys;
}


int ChiSquared::NumBin() {
	// compute the number of nonempty bins
	// if not defined, the routine will do this check
	// and creates lower and upper bin pairs in _limits
	//
	if (_nBin < 0)	// all number of bins (compressed)
		std::cerr << "WARNING - ChiSquared: binning not defined!" << std::endl;
	else if (_nBin == 0)
		std::cerr << "WARNING - ChiSquared: number of bin is 0" << std::endl;

	return _nBin;
}

void ChiSquared::DefineBinning() {
	// extract bining from spectra
	// binning is compressed, such that nonempty bins are not stored

	// unoscillated spectrum, just for counting nonempty bins
	std::map<std::string, TH1D*> spectra = BuildSpectrum();

	_nBin = 0;
	_allBin = 0;
	_limits.clear();	// for construct TH1D* / spectra
	_binpos.clear();
	_global.clear();

	for (const auto &it : spectra) {
		int b0 = std::max(it.second->FindFirstBinAbove() - 1, 0);
		int b1 = std::min(it.second->FindLastBinAbove(),
				  it.second->GetNbinsX());

		//int b1 = it.second->GetNbinsX();
		//if (it.first.find_first_of("E") != std::string::npos)
		//	b1 = it.second->FindBin(1.251);	// 1.25 GeV hard cut-off
		// store first and last bin
		// refers to TH1D
		//_limits[it.first] = b1;
		_limits[it.first] = std::pair<int, int>(b0, b1);
		// refers to _global
		_binpos[it.first] = std::pair<int, int>(_nBin,
							_nBin + b1 - b0);

		const double *bins = it.second->GetXaxis()
		/* global binning */	      ->GetXbins()->GetArray();
		//_global[it.first].assign(bins, bins + b1 + 1);
		_global[it.first].assign(bins, bins + it.second->GetNbinsX() + 1);

		// count number of total bins
		_nBin += b1 - b0;
		_allBin += it.second->GetNbinsX();

		delete it.second;
	}
	
	if (kVerbosity)
		std::cout << "Number of nonzero bins: "
			  << _nBin << " / " << _allBin << std::endl;
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

	if (!cd->Get("corr_file", matFile) && !_nScale) {
		std::cerr << "No correlation matrix file or scale error specified\n\n"
			  << "~~~ FITTING WITHOUT SYSTEMATICS ~~~\n\n";
		_nSys = 0;
		return;
	}
	
	if (matFile.empty()) {	// only scale error
		//last _nScale entries are for energy scale
		corr = Eigen::MatrixXd::Identity(_nScale, _nScale);
		_nSys = _nScale;
		return;
	}

	cd->Get("corr_name", matName);
	TFile * mf = new TFile(matFile.c_str());
	TMatrixT<double> * cmat = static_cast<TMatrixT<double>*>(mf->Get(matName.c_str()));

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	else
		sysA = std::max(0, sysA);
	if (!cd->Get("syst_last", sysB))
		sysB = cmat->GetNcols();
	else
		sysB = std::min(cmat->GetNcols(), sysB);

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
	//std::cout << "corr diag\n" << corr.diagonal() << std::endl;

	if (kVerbosity) {
		std::cout << "Importing " << sysB-sysA
			  << " entries from " << matName
			  << " ( " << matFile << " ) " << std::endl;
		std::cout << "Total number of systematics is "
			  << corr.cols() << std::endl;
	}

	_nSys = corr.cols();
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
	sysMatrix[0] = Eigen::ArrayXXd::Zero(_nBin, _nSys);
	for (int sigma = -3; sigma < 4; sigma += 2)
		sysMatrix[sigma] = Eigen::ArrayXXd::Zero(_nBin, _nSys);

	int sysA, sysB;
	if (!cd->Get("syst_first", sysA))
		sysA = 0;
	if (!cd->Get("syst_last", sysB))
		sysB = _nSys;

	int off = 0;	// offset for global bin
	for (const std::string &it : _type) {
		// it is sample type name
		std::string fileName;
		if (!cd->Get("systematic_" + it, fileName)) {
			std::cout << "ChiSquared: no systematic for " << it << " sample,"
				  << " parameters will be set to zero and not fitted\n";
			continue;
		}

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
				//for (int n = 0, j = 0; n < _limits[it]; ++n, ++j)
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
				//for (int n = 0, j = 0; n < _limits[it]; ++n, ++j)
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
		//for (int n = 0; n < _limits[is.first]; ++n, ++j) 
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
			if (kVerbosity > 2)
				std::cout << "using reconstruction file " << reco_file << std::endl;
			Reco *reco = new Reco(reco_file);

			const double *bins;
			int nBin = reco->BinsX(bins);
			f1->SetRange(bins[0], bins[nBin-1]);

			// flux contains the weight of spectrum
			std::string fname = ih + "_" + im;
			TH1D* flux = new TH1D(fname.c_str(), fname.c_str(), nBin, bins);
			flux->SetDirectory(0);
			flux->Add(f1);	//filled with ones

			if (osc) {
				if (kVerbosity > 2)
					std::cout << "Oscillating spectrum\n";
				osc->Oscillate(_fin[cm], _fout[cm], flux);
			}
			else if (kVerbosity > 2)
				std::cout << "Not oscillating\n";

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
					if (kVerbosity > 2)
						std::cout << "Projecting " << hname
							  << " from "
							  << ih << ", " << im << ", "
							  << ic << std::endl;

					if (spectra.count(hname) && spectra[hname])
						spectra[hname]->Add(py);
					else
						spectra[hname] = static_cast<TH1D*>(py->Clone());
				}
			}
			delete flux;
			delete reco;
			++cm;
		}
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
	//Eigen::VectorXd epsil = Eigen::VectorXd::Constant(_nSys, 0);
	Eigen::VectorXd best_eps = epsil;
	double best_x2 = X2(On, En, best_eps);

	int tries = 0;
	// if _nSys == 0 there is no fit, good for stats only
	while (_nSys > 0 && !FindMinimum(On, En, epsil) && tries < maxIteration) {
		++tries;
		double x2 = X2(On, En, epsil);

		if (kVerbosity > 1)
			std::cout << "~~~~~~ No convergence reached ~~~~~~~";

		double dist = (best_eps - epsil).norm();
		double step = x2 - best_x2;
		if (best_x2 > x2) {
			best_x2 = x2;
			best_eps = epsil;
			std::cout << " but new best!";
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
bool ChiSquared::FindMinimum(const Eigen::VectorXd &On,
			     const Eigen::VectorXd &En, 
			     Eigen::VectorXd &epsil)
			     //double alpha)
{
	int c = 0;	//counter
	double lambda = lm_0;

	Eigen::VectorXd delta = Eigen::VectorXd::Constant(_nSys, 1);

	double x2 = X2(On, En, epsil);
	double diff = 1;
	//double minj = 1;
	//while (diff > err && minj > NumSys() * err && minj > 1e-6 &&
	//       lambda < 1e6 && c < maxIteration)
	//while (delta.norm() > NumSys() * err && c < maxIteration)

	if (kVerbosity > 2)
		std::cout << "Minimising fit from x2: " << x2 << std::endl;

	while (std::abs(diff) / DOF() > fitErr
	      && delta.norm() / _nSys > fitErr) {
		++c;	//counter

		// build hessian and gradient/jacobian together
		// to save computation time
		Eigen::VectorXd jac(_nSys);
		Eigen::MatrixXd hes(_nSys, _nSys);
		JacobianHessian(jac, hes, On, En, epsil);

		//std::cout << "jac\n" << jac.transpose() << "\n"
		//	  << "hes\n" << hes.transpose() << std::endl;
		//add diagonal to hesT hes
		double maxd = hes.diagonal().maxCoeff();
		hes.diagonal() += Eigen::VectorXd::Constant(_nSys, maxd * lambda);
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
		else if (std::abs(delta.norm()) > 0)	//next x2 is worse
			lambda *= lm_up;	//nothing changes but lambda
						//if step is nonnull

		//std::cout << "\ndelta\n" << delta << std::endl;
		//if (false)
		if (kVerbosity > 2) {
			std::cout << c << " -> l " << lambda
				  << ",\tstep: " << std::abs(delta.norm() / NumSys())
				  << ", X2: " << x2
				  << " ( " << std::abs(diff / DOF())
				  << " ) " << std::endl;
		}
	}


	return (lambda <= 1./lm_0);	//convergence was reached
					// else	no convergence, change starting point
}


Eigen::VectorXd ChiSquared::Jacobian(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En, 
				      const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(_nSys);
	Eigen::MatrixXd hes(_nSys, _nSys);
	JacobianHessian(jac, hes, On, En, epsil);

	return jac;
}


Eigen::MatrixXd ChiSquared::Hessian(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En, 
				     const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(_nSys);
	Eigen::MatrixXd hes(_nSys, _nSys);
	JacobianHessian(jac, hes, On, En, epsil);

	return hes;
}

/*
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
		std::string it = TypeFromBin(n);
		int t = _nSys - _nScale + _type_scale[it];

		double en = Scale(Ep, sk, n, it);	// en can be empty!
		if (en < 1e-9)			//bad stuff, avoid it
			continue;

		double one_oe = 1 - On(n) / en;
		double en_jac = ScaleJac(Ep, sk, n, it);

		// hessian of scale error first
		hes(t, t) += On(n) * pow(en_jac / en, 2);
		//std::cout << "adding " << n << "\t"
			  //<< On(n) * pow(en_jac / en, 2)
			  //<< " += " << hes(t, t) << std::endl;;

		if (std::abs(one_oe) > 1e-9) {
			jac(t)    += one_oe * en_jac;
			hes(t, t) += one_oe * ScaleHes(Ep, sk, n, it);
		}

		//do the rest, including mixed term with SK syst
		for (int k = 0; k < _nSys - _nScale; ++k) {
			const Eigen::VectorXd &gk = gam.col(k);
			double gn_jac = Scale(gk, sk, n, it);

			//mixed term with SK error
			hes(k, t) += en_jac * gn_jac * On(n) / en / en;

			//diagonal term
			hes(k, k) += On(n) * pow(gn_jac / en, 2);

			if (std::abs(one_oe) < 1e-9) {
				jac(k)    += one_oe * gn_jac;
				hes(k, t) += one_oe * ScaleJac(gk, sk, n, it);
			}

			for (int j = k + 1; j < _nSys - _nScale; ++j) {
				const Eigen::VectorXd &gj = gam.col(j);

				hes(k, j) += gn_jac * Scale(gj, sk, n, it)
					   * On(n) / en / en;

				if (std::abs(one_oe) < 1e-9)
					continue;

				const Eigen::VectorXd gh = GammaJac(gk, epsil, j);
				hes(k, j) += one_oe * Scale(gh, sk, n, it);
			}
		}
	}

	for (int k = 0; k < _nSys - _nScale; ++k)
		for (int j = k+1; j < _nSys - _nScale; ++j)
			hes(j, k) = hes(k, j);
}



void ChiSquared::JacobianHessian3(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				    const Eigen::VectorXd &On,
				    const Eigen::VectorXd &En, 
				    const Eigen::VectorXd &epsil)
{
	//corr is inverse of correlation matrix
	jac = corr * epsil;	//gradient/jacobian
	//hes = corr;		//hessian
	hes.setZero();		//hessian

	// just scale energy bit
	const Eigen::VectorXd sk  = epsil.tail(_nScale);

	// event distribution with systematics applied
	const Eigen::VectorXd Ep  = Gamma(En, epsil);

	// matrix contains derivative terms for each systematics and for each bin
	const Eigen::MatrixXd Fp = gp(epsil);
	const Eigen::VectorXd Uno = Eigen::VectorXd::Constant(_nBin, 1);

	// loop over samples
	for (const std::string &it : _type) {

		// scale error value
		double scale_err = _scale[it];
		// energy shift
		double shift = 1 + sk(_type_scale[it]) * scale_err;
		// offset between global bin and energy bin
		int off = _binpos[it].first - _limits[it].first;
		// absolute systematic error for this scale parameter
		int t = _nSys - _nScale + _type_scale[it];

		// alias to binning
		std::vector<double> &global = _global[it];

		// loop over bins of this sample
	for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {

		// scale systematic is related to bin number
		// scaling bin start
		int m0 = StartingBin(it, shift, n);
		int m1 = global.size() - 1;
		//int m1 = EndingBin(it, shift, n);
		// unscaled/original bins
		double b0_n = global[n - off];
		double b1_n = global[n - off + 1];

		//std::cout << "bin " << n << " in shifted " << m0 << ", " << m1
		//	  << "\n\t" << global[m0] * shift
		//	  << " - " << global[m1-1] * shift 
		//	  << "\t" << b0_n << " - " << b1_n << std::endl;


		// shifted events in this bin, it can never be empty
		double en = QuickScale(&ChiSquared::Nor, Ep,
				global, b0_n, b1_n,
				scale_err, shift, m0, m1, off);

		double one_oe = 1 - On(n) / en;
		double en_jac = QuickScale(&ChiSquared::Jac, Ep,
				global, b0_n, b1_n,
				scale_err, shift, m0, m1, off);

		// jacobian
		jac(t)    += one_oe * en_jac;

		// hessian of scale error first
		hes(t, t) += On(n) * pow(en_jac / en, 2)
			+ one_oe * QuickScale(&ChiSquared::Hes, Ep,
					global, b0_n, b1_n,
					scale_err, shift, m0, m1, off);

		//do the rest, including mixed term with SK syst
		for (int k = 0; k < _nSys - _nScale; ++k) {
			// scaled jacobian in this bin
			Eigen::VectorXd kk = Ep.cwiseProduct(Fp.col(k));
			double gn_jac = QuickScale(&ChiSquared::Nor, kk,
					global, b0_n, b1_n,
					scale_err, shift, m0, m1, off);

			// jacobian
			jac(k)    += one_oe * gn_jac;

			//diagonal term
			hes(k, k) += On(n) * pow(gn_jac / en, 2);

			//mixed term with SK error
			hes(k, t) += en_jac * gn_jac * On(n) / en / en
				+ one_oe * QuickScale(&ChiSquared::Jac, kk,
						global, b0_n, b1_n,
						scale_err, shift, m0, m1, off);

			// mixed terms with non energy scale parameters
			for (int j = k + 1; j < _nSys - _nScale; ++j) {

				//std::cout << "n " << n << "\tk " << k
				//	  << "\tj " << j << "\n";
				Eigen::VectorXd jj = Ep.cwiseProduct(Fp.col(j));
				Eigen::VectorXd kj = kk.cwiseProduct(Fp.col(j));
				hes(k, j) += gn_jac * On(n) / en / en
					* QuickScale(&ChiSquared::Nor, jj,
							global, b0_n, b1_n,
							scale_err, shift, m0, m1, off)
					+ one_oe * QuickScale(&ChiSquared::Nor, kj,
							global, b0_n, b1_n,
							scale_err, shift, m0, m1, off);
			}
		}
		}
	}

	
	// make matrix symmetric
	for (int k = 0; k < _nSys - _nScale; ++k)
		for (int j = k+1; j < _nSys - _nScale; ++j)
			hes(j, k) = hes(k, j);
}
*/

void ChiSquared::JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				 const Eigen::VectorXd &On,
				 const Eigen::VectorXd &En, 
				 const Eigen::VectorXd &epsil)
{
	//corr is inverse of correlation matrix
	jac = corr * epsil;	//gradient/jacobian
	hes = corr;		//hessian
	//hes.setZero();		//hessian

	// event distribution with systematics applied
	const Eigen::ArrayXd Ep = Gamma(En, epsil);
	// matrix contains derivative terms for each systematics and for each bin
	const Eigen::ArrayXXd Fp = one_Fp(epsil);

	// loop over samples
	for (const std::string &it: _type) {
		int t = _nScale ? _nSys - _nScale + _type_scale[it] : 0;
		double skerr = t < epsil.size() ? epsil(t) : 0;

		std::vector<std::pair<int, int> > slices = AllSlices(it, skerr);
		std::vector<Eigen::ArrayXd> scales = AllScale(&ChiSquared::Nor,
								it, skerr);
		std::vector<Eigen::ArrayXd> jacobs = AllScale(&ChiSquared::Jac,
								it, skerr);
		std::vector<Eigen::ArrayXd> hessis = AllScale(&ChiSquared::Hes,
								it, skerr);

	for (int n = 0; n < _binpos[it].second - _binpos[it].first; ++n) {
	//for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {

		int m0 = slices[n].first, dm = slices[n].second - m0;
		const Eigen::ArrayXd ev = Ep.segment(m0, dm);

		if (kVerbosity > 3)
			std::cout << "type " << it << ", bin "
				  << n + _binpos[it].first
				  << " scales are between "
				  << m0 << " and " << m0 + dm << std::endl;

		double on = On(n + _binpos[it].first);
		double en = (scales[n] * ev).sum();

		double one_oe = 1 - on / en;

		double en_jac = (jacobs[n] * ev).sum();

		// jacobian
		if (_nScale) {
			jac(t)    += one_oe * en_jac;

			// hessian of scale error first
			//std::cout << "large " << on << ", " << en
			//	  << ", " << en_jac << ", " << one_oe
			//	  << ";\n\t" << ev.transpose()
			//	  << ";\n\t" << scales[n].transpose()
			//	  << ";\n\t" << jacobs[n].transpose()
			//	  << ";\n\t" << hessis[n].transpose()
			//	  << " = " << (hessis[n] * ev).sum() << std::endl;
			//std::cout << "->terms " << on * pow(en_jac / en, 2)
			//	  << "\t" << one_oe * (hessis[n] * ev).sum() << std::endl;
			hes(t, t) += on * pow(en_jac / en, 2)
				+ one_oe * (hessis[n] * ev).sum();
		}

		//do the rest, including mixed term with SK syst
		for (int k = 0; k < _nSys - _nScale; ++k) {
			// scaled jacobian in this bin
			Eigen::ArrayXd kk = Fp.col(k).segment(m0, dm);

			double gn_jac = (scales[n] * ev * kk).sum();

			// jacobian
			jac(k)    += one_oe * gn_jac;

			//diagonal term
			hes(k, k) += on * pow(gn_jac / en, 2);

			if (_nScale)
				//mixed term with SK error
				hes(k, t) += en_jac * gn_jac * on / en / en
					+ one_oe * (jacobs[n] * kk).sum();

			// mixed terms with non energy scale parameters
			for (int j = k + 1; j < _nSys - _nScale; ++j) {

				if (kVerbosity > 5)
					std::cout << "n " << n << "\tk " << k
						  << "\tj " << j << "\n";
				Eigen::ArrayXd jj = Fp.col(j).segment(m0, dm);

				hes(k, j) += gn_jac * on / en / en
					  * (scales[n] * ev * jj).sum()
					+ one_oe
					* (scales[n] * ev * jj * kk).sum();
			}
		}
	}
	//std::cout << "\tend type\n";
	}

	// make matrix symmetric
	hes.triangularView<Eigen::Lower>() = hes.transpose();
}


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
double ChiSquared::ObsX2(const Eigen::VectorXd &On,
			 const Eigen::VectorXd &En,
			 const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	return ObsX2n(On, En, epsil).sum();
}


// return statistics X2 as a vector
// so X2 contribution to each bin
Eigen::ArrayXd ChiSquared::ObsX2n(const Eigen::VectorXd &On,
				  const Eigen::VectorXd &En,
				  const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	Eigen::ArrayXd gam = Gamma(En, epsil);

	Eigen::ArrayXd chi2(_nBin);
	for (const std::string &it : _type) {
		int t = _nScale ? _nSys - _nScale + _type_scale[it] : 0;
		double skerr = t < epsil.size() ? epsil(t) : 0;
		for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {
			double en = Scale(gam, skerr, n, it);
			chi2(n) = 2 * en - 2 * On(n) * (1 + log(en / On(n)));
		}
	}

	return chi2;
}


// return systematic X2
double ChiSquared::SysX2(const Eigen::VectorXd &epsil) {
	return epsil.transpose() * corr * epsil;
}


// this return the spectrum modified with the systematics
Eigen::ArrayXd ChiSquared::Gamma(const Eigen::VectorXd &En,
				 const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXd gam = En.matrix();
	for (int k = 0; k < _nSys - _nScale; ++k)
		gam *= one_Fk(epsil(k), k);

	return gam;
}


/*
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
*/

// this function return a modifier for the bin relative to jacobians of the gamma function
// it returns a correct factor only if a valid bin is passed,
// otherwise it returns 1 as backward-compatibility
//double ChiSquared::GammaHes(int n, int k, double ek, int j, double ej)
//{
//	if (k < 0)
//		return 1;
//	else if (j < 0)
//		return gp(k, n, ek);
//	else
//		return gp(k, n, ek) * gp(j, n, ej);
//}


// this is (1 + F)
Eigen::ArrayXXd ChiSquared::one_F(const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXXd f(_nBin, _nSys);
	
	for (int k = 0; k < epsil.size(); ++k)
		f.col(k) = one_Fk(epsil(k), k);

	return f;
}

// k-th entry of the previous function
Eigen::ArrayXd ChiSquared::one_Fk(double err, int k)
{
	double dl, du;

	if (err < -1) {
		dl = -3; du = -1;
	}
	else if (err >= -1 && err < 0) {
		dl = -1; du = 0;
	}
	else if (err >= 0 && err < 1) {
		dl = 0; du = 1;
	}
	else if (err >= 1) {
		dl = 1; du = 3;
	}

	const Eigen::ArrayXd &sl = sysMatrix[dl].col(k);
	Eigen::ArrayXd su = sysMatrix[du].col(k);

	su = (su - sl) / (du - dl);
	return Eigen::ArrayXd::Constant(_nBin, 1) + su * (err - dl) + sl;
}

// this is Fp / (1 + F) which is derivative wrt epsilon
Eigen::ArrayXXd ChiSquared::one_Fp(const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXXd fp(_nBin, _nSys);
	
	for (int k = 0; k < epsil.size(); ++k)
		fp.col(k) = one_Fpk(epsil(k), k);

	return fp;
}

// k-th entry of the previous function
Eigen::ArrayXd ChiSquared::one_Fpk(double err, int k)
{
	double dl, du;

	if (err < -1) {
		dl = -3; du = -1;
	}
	else if (err >= -1 && err < 0) {
		dl = -1; du = 0;
	}
	else if (err >= 0 && err < 1) {
		dl = 0; du = 1;
	}
	else if (err >= 1) {
		dl = 1; du = 3;
	}

	Eigen::ArrayXd sl = sysMatrix[dl].col(k);
	Eigen::ArrayXd su = sysMatrix[du].col(k);

	su = (su - sl) / (du - dl);
	sl = Eigen::ArrayXd::Constant(_nBin, 1) + su * (err - dl) + sl;

	return su / sl;
}

/*
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
*/



// Apply energy scale dilation
// FactorFn is a function which changes if calculating energy scale, jacobian or hessian
// The dilation does not depend on En, therefore 
// En can be anything
//double ChiSquared::Scale(FactorFn factor, const Eigen::VectorXd &En,
//			const Eigen::VectorXd &sigma, int t)

std::string ChiSquared::TypeFromBin(int n) {
	// check if is the same type
	if (n >= _binpos[_tlim].first && n < _binpos[_tlim].second)
		return _tlim;

	// return sample type given bin
	for (const auto &it : _binpos)
		if (it.first != _tlim
		    && n >= it.second.first
		    && n < it.second.second) {
			_tlim = it.first;
			return _tlim;
		}

	return "";
}

int ChiSquared::StartingBin(std::string it, double shift, int n)
{
	n += _limits[it].first - _binpos[it].first;
	auto im = std::lower_bound(_global[it].begin(),
				   _global[it].end(),
				   _global[it][n] / shift);
	return std::max(int(std::distance(_global[it].begin(), im)) - 1,
			_limits[it].first);
}

int ChiSquared::EndingBin(std::string it, double shift, int n)
{
	n += _limits[it].first - _binpos[it].first;
	auto im = std::upper_bound(_global[it].begin(),
				   _global[it].end(),
				   _global[it][n+1] / shift);
	return std::min(int(std::distance(_global[it].begin(), im)),
			_limits[it].second);
}

std::vector<std::pair<int, int> > ChiSquared::AllSlices(std::string it,
							double skerr) {

	std::vector<std::pair<int, int> > allslices;

	// energy shift
	double shift = 1 + skerr * _scale[it];
	// offset between global bin and energy bin
	int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// loop over bins of this sample
	for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {
		//std::cout << "sbin " << n << " is "
		//	<< global[n - off] << " - " << global[n+1-off] << std::endl;

		// scale systematic is related to bin number
		// scaling bin start

		//std::cout << "\tshift is " << shift * global[m0] << " (" << m0+off
		//	  << ")  - " << shift * global[m1] << " (" << m1+off << ")\n";

		if (_nScale) {
			int m0 = StartingBin(it, shift, n);
			int m1 = EndingBin(it, shift, n);
			allslices.push_back(std::make_pair(m0 + off, m1 + off));
		}
		else
			allslices.push_back(std::make_pair(n, n+1));
	}

	return allslices;
}


std::vector<Eigen::ArrayXd> ChiSquared::AllScale(FactorFn factor,
						 std::string it, double skerr)
{
	std::vector<Eigen::ArrayXd> allfacts;

	if (!_nScale) {	// no scaling needed
		for (int n = _binpos[it].first; n < _binpos[it].second; ++n)
			allfacts.push_back(Eigen::ArrayXd::Constant(1, 1));
		return allfacts;
	}

	// scale error value
	double scale_err = _scale[it];
	// energy shift
	double shift = 1 + skerr * scale_err;
	// offset between global bin and energy bin
	int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// alias to binning
	std::vector<double> &global = _global[it];

	// loop over bins of this sample
	for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {

		// scale systematic is related to bin number
		// scaling bin start
		int m0 = StartingBin(it, shift, n);
		int m1 = EndingBin(it, shift, n);

		// unscaled/original bins
		double b0_n = global[n - off];
		double b1_n = global[n - off + 1];

		std::vector<double> thisfact;
		for (int m = m0; m < m1; ++m) {

			// scaled bins
			double b0_m = global[m];
			double b1_m = std::min(global[m + 1], 30 / shift);
			//std::cout << "contains " << b0_m << " and " << b1_m << std::endl;
			//continue/break computation because there is no overlap
			if (shift * b0_m > b1_n) {
				//std::cout << "stopping at " << m
				//	  << " instead of " << m1-1 << std::endl;
				break;
			}

			double f = (b1_n - b0_n + shift * (b1_m - b0_m)
					- std::abs(b0_n - shift * b0_m)
					- std::abs(b1_n - shift * b1_m)) / 2.;

			int s0 = b0_n - shift * b0_m < 0 ? -1 : 1;
			int s1 = b1_n - shift * b1_m < 0 ? -1 : 1;
			double ss = f > 0 ? 1 : 0.5;	//continuity factor
			double fd = ss * (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;

			std::vector<double> terms = {scale_err, shift,
						     f, fd, b1_m - b0_m};
			thisfact.push_back((this->*factor)(terms));
		}

		Eigen::ArrayXd arr = Eigen::Map<Eigen::ArrayXd>
					(thisfact.data(), thisfact.size(), 1);
		allfacts.push_back(arr);
	}

	return allfacts;
}

/*
double ChiSquared::QuickScale(FactorFn factor,
			 const Eigen::VectorXd &En,
			 //const Eigen::VectorXd &gh,
			 const std::vector<double> &global,
			 double b0_n, double b1_n,
			 double scale_err, double shift,
			 int m0, int m1, int off)
			 //int k, double ek, int j, double ej)
{
	double ret = 0;
	// Loop over unscaled/original edges only nonempty bins
	for (int m = m0; m < m1; ++m) {
		// scaled bins
		double b0_m = global[m];
		double b1_m = std::min(global[m + 1], 30 / shift);
		//continue/break computation because there is no overlap
		if (shift * b0_m > b1_n)
			break;

		double f = (b1_n - b0_n + shift * (b1_m - b0_m)
			  - std::abs(b0_n - shift * b0_m)
			  - std::abs(b1_n - shift * b1_m)) / 2.;

		int s0 = b0_n - shift * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - shift * b1_m < 0 ? -1 : 1;
		double fd = (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;
		//double ss = f > 0 ? 1 : 0.5;	//continuity factor

		std::vector<double> terms = {scale_err, shift,
					     f, fd, b1_m - b0_m};
		ret += En(m + off) * (this->*factor)(terms);
	}

	return ret;
}
*/

double ChiSquared::Scale(FactorFn factor,
			 const Eigen::ArrayXd &En,
			 double skerr, int n, std::string it, int m0)
{
	if (!_nScale) 	// no scaling
		return En(n);

	if (it.empty())
		it = TypeFromBin(n);

	double scale_err = _scale[it];	// this is Error value
	double shift = 1 + skerr * scale_err;
	//
	// binary search for starting bin
	// return starting point for following loop

	if (m0 < 0)
		m0 = StartingBin(it, shift, n);

	int off = _binpos[it].first - _limits[it].first;

	// unscaled/original bins
	double b0_n = _global[it][n - off];
	double b1_n = _global[it][n - off + 1];

	//std::cout << "type " << it << " from " << _binpos[it].first << " to " << _binpos[it].second << std::endl;
	//std::cout << "bin " << n << " : " << b0_n << " - " << b1_n << std::endl;
	//std::cout << "starting from " << m0 << ": "
		  //<< shift * _global[it][m0] << " and " << shift * _global[it][m0 + 1] << "\n";
	//std::cout << "offset " << off << " from " << m0 + off << std::endl;

	double ret = 0;
	// Loop over unscaled/original edges only nonempty bins
	for (int m = m0; m < _global[it].size()-1; ++m) {
		// scaled bins
		double b0_m = _global[it][m];
		double b1_m = std::min(_global[it][m + 1], 30 / shift);
		//continue/break computation because there is no overlap
		if (shift * b0_m > b1_n)
			break;

		double f = (b1_n - b0_n + shift * (b1_m - b0_m)
			  - std::abs(b0_n - shift * b0_m)
			  - std::abs(b1_n - shift * b1_m)) / 2.;

		//if (f < 0)	// should never be negative
		//	continue;


		int s0 = b0_n - shift * b0_m < 0 ? -1 : 1;
		int s1 = b1_n - shift * b1_m < 0 ? -1 : 1;
		double ss = f > 0 ? 1 : 0.5;	//continuity factor
		double fd = ss * (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;

		std::vector<double> terms = {scale_err, shift, f, fd, b1_m - b0_m};
		ret += En(m + off) * (this->*factor)(terms);
		//std::cout << "  " << m << " : "
		//	<< b0_m * shift << " - "
		//	<< b1_m * shift << "\t" << En(m + off)
		//	  << "\t" << (this->*factor)(terms) << std::endl;

	}
	//std::cout << "returning " << ret << std::endl << std::endl;

	return ret;
}


// call Scale Energy function
double ChiSquared::Scale(const Eigen::VectorXd &En,
			 double skerr,
			 int n, std::string it, int m0)
{
	const Eigen::ArrayXd &arrEn = En.array();
	return Scale(&ChiSquared::Nor, arrEn, skerr, n, it, m0);
}

double ChiSquared::Scale(const Eigen::ArrayXd &En,
			 double skerr,
			 int n, std::string it, int m0)
			 //int k, double ek, int j, double ej)
{
	return Scale(&ChiSquared::Nor, En, skerr, n, it, m0); //, k, ek, j, ej);
}

/*
// call jacobian of Scale Energy function
double ChiSquared::ScaleJac(const Eigen::VectorXd &En,
			 double skerr,
			 int n, std::string it, int m0)
{
	Eigen::ArrayXd arrEn = En.array();
	return Scale(&ChiSquared::Nor, arrEn, skerr, n, it, m0);
}

double ChiSquared::ScaleJac(const Eigen::VectorXd &En,
			    double skerr,
			    int n, std::string it, int m0)
			    //int k, double ek, int j, double ej)
{
	return Scale(&ChiSquared::Jac, En, skerr, n, it, m0); //, k, ek, j, ej);
}


// call hessian of Scale Energy function
double ChiSquared::ScaleHes(const Eigen::VectorXd &En,
			    double skerr,
			    int n, std::string it, int m0)
			    //int k, double ek, int j, double ej)
{
	return Scale(&ChiSquared::Hes, En, skerr, n, it, m0); //, k, ek, j, ej);
}
*/


double ChiSquared::Nor(const std::vector<double> &term)
{
	//return f / shift / db;
	return term[2] / term[1] / term[4];
}

double ChiSquared::Jac(const std::vector<double> &term)
{
	//return scale_err / shift / db * (fd - f / shift)
	return term[0] / term[1] / term[4] * (term[3] - term[2] / term[1]);
}

double ChiSquared::Hes(const std::vector<double> &term)
{
	//return pow(scale_err / shift, 2) / db * (f / shift - fd)
	return 2 * pow(term[0] / term[1], 2) / term[4]
		 * (term[2] / term[1] - term[3]);
}
