#include "BeamSample.h"

//ctor
BeamSample::BeamSample(const std::string &card, std::string process) :
	Sample(card)
{
	CardDealer cd(card);
	Init(cd, process);
}

BeamSample::BeamSample(const CardDealer &cd, std::string process) :
	Sample(cd)
{
	Init(cd, process);
}


BeamSample::BeamSample(CardDealer *cd, std::string process) :
	Sample(cd)
{
	Init(*cd, process);
}

void BeamSample::Init(const CardDealer &cd, std::string process)
{
	//_mode = {"nuE0_nuE0", "nuM0_nuM0", "nuM0_nuE0",
	//	 "nuEB_nuEB", "nuMB_nuMB", "nuMB_nuEB"};
	//_chan = {"E_CCQE", "M_CCQE", "E_CCnQE", "M_CCnQE", "E_NC", "M_NC"};
	//_horn = {"FHC", "RHC"};
	//std::map<std::string, std::pair<flavour, flavour> >
	_oscf = { {"nuE0_nuE0", {Nu::E_, Nu::E_}},
		  {"nuM0_nuM0", {Nu::M_, Nu::M_}},
		  {"nuM0_nuE0", {Nu::M_, Nu::E_}},
		  {"nuEB_nuEB", {Nu::Eb, Nu::Eb}},
		  {"nuMB_nuMB", {Nu::Mb, Nu::Mb}},
		  {"nuMB_nuEB", {Nu::Mb, Nu::Eb}} };
	// _fin = {Nu::E_, Nu::M_, Nu::M_, Nu::Eb, Nu::Mb, Nu::Mb};
	//_fout = {Nu::E_, Nu::M_, Nu::E_, Nu::Eb, Nu::Mb, Nu::Eb};

	// if sample beam is not defined..
	if (!cd.Get("sample_beam", _type))
		// default value
		_type = {"E_FHC", "E_RHC", "M_FHC", "M_RHC"};

	// get scale error which is now defined manually
	if (!cd.Get("scale_error_", _scale)) {
		// check first for a map, if not look for a scalar
		double err;
		_scale.clear();
		if (cd.Get("scale_error", err)) {
			_nScale = 1;
			for (const std::string &it : _type) {
				_scale[it] = err;
				_type_scale[it] = 0;
			}
		}
		else
			_nScale = 0;
	}
	else {	// there should be scale_E and scale_M
		// map scale contains "E" and "M"	
		std::map<std::string, double> err;
		for (const std::string &it : _type) {
			//E_FHC or E_RHC
			if (it.find_first_of('E') != std::string::npos) {
				if (_scale.count("E"))
					err[it] = _scale["E"];
				else {
					if (kVerbosity)
						std::cout << "No scale_E in card"
							  << cd.CardName() << std::endl;
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
							  << cd.CardName() << std::endl;
					// default to
					err[it] = 0.024;
				}
				_type_scale[it] = 1;
			}
		}
		_nScale = _scale.size();	//of size 1, 2, or 4
		_scale = err;			// reassign new map
	}


	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "BeamSample: types: ";
		for (const auto &it : _type)
			std::cout << "\t" << it;
		std::cout << std::endl;
		std::cout << "BeamSample: number of scale systematics " << _nScale << std::endl;
		for (const auto &it : _scale)
			std::cout << "\t" << it.first << " -> " << it.second << std::endl;
	}

	std::string profile;
	if (!cd.Get("density_profile", profile) && kVerbosity)
		std::cout << "BeamSample: density profile not set in card, please set it yourself" << std::endl;
	else
		_lens_dens = Oscillator::GetMatterProfile(profile);



	if (process.empty())
		process = "RBS";
	else
		std::transform(process.begin(), process.end(), process.begin(),
				[](unsigned char c){ return std::toupper(c); });

	if (process.find('R') != std::string::npos)
		LoadReconstruction(cd);
	if (process.find('B') != std::string::npos)
		DefineBinning();	// from Sample.h
	if (process.find('S') != std::string::npos)
		LoadSystematics(cd);
}


void BeamSample::LoadReconstruction(const CardDealer &cd)
{
	// open it if already exists
	std::vector<std::string> reco_files;
	if (!cd.Get("reco_input", reco_files))
		throw std::invalid_argument("BeamSample: no reconstruction files in card,"
					    "very bad!");
	for (const std::string &reco : reco_files)
		LoadReconstruction(reco);

//	for (const std::string &ih : _horn)	//FHC, RHC
//		for (const std::string &im : _mode)	//nuE->nuE, nuM->nuE, nuM->nuM
//			LoadReconstruction(im + "_" + ih);
			//LoadReconstruction("beam_" + im + "_" + ih);

	//cmd = "rm .reconstruction_files";
	//system(cmd.c_str());
}

void BeamSample::LoadReconstruction(std::string reco_file)
{
//	std::string reco_file, type = "reco_" + channel;
//	if (!cd->Get(type, reco_file)) {
//		if (kVerbosity > 1)
//			std::cout << "WARNING - BeamSample: reconstruction file for channel \""
//				  << reco_file << "\" is missing" << std::endl;
//		return;
//	}
//

	if (reco_file.find("nu") == std::string::npos) {
		std::cout << "BeamSample: file " << reco_file << " does not have "
			<< "channel type (nu.._nu.._.HC format) in its name" << std::endl;
		return;
	}	

	std::string channel = reco_file.substr(reco_file.find("nu"));
	if (channel.find(".card") != std::string::npos)
		channel.erase(channel.find(".card"));

	if (kVerbosity > 2)
		std::cout << "BeamSample: using reconstruction file " << reco_file
			  << " for channel " << channel << std::endl;

	CardDealer cd(reco_file);

	std::string reco_path;
	if (!cd.Get("reco_path", reco_path)) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - BeamSample: no reco_path specified in \""
				  << reco_file << "\"" << std::endl;
		return;
	}

	double weight;
	if (!cd.Get("scale", weight))
		weight = 1.0;

	//std::cout << "opening " << fileReco << std::endl;
	TFile* inFile = new TFile(reco_path.c_str(), "READ");
	if (inFile->IsZombie()) {
		if (kVerbosity > 0)
			std::cerr << "WARNING - BeamSample: file " << reco_path
				  << " does not exist" << std::endl;
		return;
	}

	std::map<std::string, std::string> mSD;
	//look if there is a particular event
	if (cd.Get("ring_", mSD))
		for (const auto &is : mSD) {
			std::string type = is.first.substr(0, is.first.find_first_of('_'))
					 + channel.substr(channel.find_last_of('_'));
			if (_type.find(type) == _type.end())
				continue;


			// TH2D is made of (xbins + 2) ( ybins + 2)
			// y0 : x0 x1 x2 x3 x4 ...
			// y1 : x0 x1 x2 x3 x4 ...
			// y2 : x0 x1 x2 x3 x4 ...
			// y3 : x0 x1 x2 x3 x4 ...
			// y4 : x0 x1 x2 x3 x4 ...
			// .....
			// as there are overflow and underflow bins
			// correct array starts from 1 to ybins + 1
			if (!inFile->Get(is.second.c_str()))
				continue;

			TH2D *h2 = static_cast<TH2D*>
				(inFile->Get(is.second.c_str()));
		
			int xs = h2->GetNbinsX();// E_true and n_cols
			int ys = h2->GetNbinsY();// E_reco and n_rows
			// catching xs+2 X y2 section
			Eigen::MatrixXd rm = Eigen::Map<Eigen::MatrixXd>
				(h2->GetArray() + xs + 2, xs + 2, ys);
			//std::cout << "raw\n";
			//for (int g = 0; g < 801; ++g)
			//	std::cout << *(h2->GetArray()+g) << "\t";
			//std::cout << std::endl;

			// remove first and last columns
			rm.transposeInPlace();
			rm.leftCols(xs) = rm.middleCols(1, xs);
			rm.conservativeResize(ys, xs);

			rm *= weight;
			// style is E_CCQE_nuM0_nuM0_RHC 
			_reco[is.first + "_" + channel] = rm;

			const double *bx = h2->GetXaxis()->GetXbins()->GetArray();
			//const double *by = h2->GetYaxis()->GetXbins()->GetArray();
			if (!_global.count(type))
				_global[type].assign(bx, bx + xs + 1);
			//_binX[is.first].assign(bx, bx + xs + 1);
			//_binY[is.first].assign(by, by + ys + 1);
		}

	inFile->Close();
}


std::map<std::string, Eigen::VectorXd>
	BeamSample::BuildSamples(std::shared_ptr<Oscillator> osc)
{
	if (osc)
		osc->SetMatterProfile(_lens_dens);

	std::map<std::string, Eigen::VectorXd> samples;
	for (const auto &ir : _reco) {

		if (kVerbosity > 4)
			std::cout << "BeamSample: using " << ir.first << std::endl;
		if (ir.first.find("NC") != std::string::npos &&
		     ( ir.first.find("nuM0_nuE0") != std::string::npos
		    || ir.first.find("nuMB_nuEB") != std::string::npos ) ) {
			std::cout << "skip " << ir.first << std::endl;
			continue;	// this shouldn't exist, but just in case...
		}

		std::string hname = ir.first;
		size_t len = hname.find_last_of('_') - hname.find_first_of('_');

		// this hould be 'nuE0_nuE0' like
		std::string cname = hname.substr(hname.find("nu"), 9);
		// this hould be 'E_FHC' like
		hname.erase(hname.find_first_of('_'), len);

		Eigen::VectorXd probs = Eigen::VectorXd::Ones(ir.second.cols());
		if (osc && ir.first.find("NC") == std::string::npos) {
			// not a NC channel
			if (kVerbosity > 4)
				std::cout << "Oscillating spectrum for " << ir.first
					  << " (" << hname << ") with "
					  << cname << "\n";
			const auto &chan = _oscf[cname];
			const auto &bins = _global[hname]; //type?

			Eigen::VectorXd vb(bins.size() - 1);
			for (int i = 0; i < vb.size(); ++i) {
				probs(i) = osc->Probability(chan.first, chan.second,
						(bins[i] + bins[i+1]) / 2.);
			}

			//probs = osc->Oscillate(chan.first, chan.second, _global[ir.first]);
		}
		else if (kVerbosity > 4)
			std::cout << "Not oscillating spectrum for " << ir.first
				<< " (" << hname << ") with "
				<< cname << "\n";


		if (samples.count(hname))
			samples[hname] += ir.second * probs;
		else
			samples[hname] = ir.second * probs;
	}

	return samples;
}



/*
 * This routine finds the systematic files and load them into a collection of Matrices.
 * The collection is for spline implementation: it is a mapping between sigma value (int) 
 * and the systematic matrix at that sigma value of the error ( NumBin * NumSys ).
 * If no systematic file is specified in the card, then the matrices are all empty 
 * and it should be like fitting for the stats only, i.e. no fit at all.
 */
void BeamSample::LoadSystematics(const CardDealer &cd)
{
	// but first, let me take the correlation matrix!

	std::string mat_file, mat_name;
	if (!cd.Get("corr_file", mat_file) && !_nScale) {
		std::cerr << "No correlation matrix file or scale error specified\n\n"
			  << "~~~ FITTING WITHOUT SYSTEMATICS ~~~\n\n";
		_nSys = 0;
		return;
	}
	
	if (mat_file.empty()) {	// only scale error
		//last _nScale entries are for energy scale
		_corr = Eigen::MatrixXd::Identity(_nScale, _nScale);
		_nSys = _nScale;
		return;
	}

	if (!cd.Get("corr_name", mat_name))
		mat_name = "correlation";

	TFile * mf = new TFile(mat_file.c_str());
	if (mf->IsZombie())
		throw std::invalid_argument("BeamSample: cannot open file for beam systematics,"
			    "cannot determine how many errors");

	TMatrixT<double> * cmat = static_cast<TMatrixT<double>*>(mf->Get(mat_name.c_str()));

	int sysA, sysB;
	if (!cd.Get("syst_first", sysA))
		sysA = 0;
	else
		sysA = std::max(0, sysA);
	if (!cd.Get("syst_last", sysB))
		sysB = cmat->GetNcols();
	else
		sysB = std::min(cmat->GetNcols(), sysB);

	if (kVerbosity)
		std::cout << "BeamSample: Limiting systematics between "
			  << sysA << " and " << sysB << std::endl;
	//last _nScale entries are for energy scale
	_corr = Eigen::MatrixXd::Zero(_nScale + sysB - sysA, _nScale + sysB - sysA);
	_corr.bottomRightCorner(_nScale, _nScale) = Eigen::MatrixXd::Identity(_nScale, _nScale);

	for (int r = sysA; r < sysB; ++r)
		for (int c = sysA; c < sysB; ++c)
			_corr(r + - sysA, c + - sysA) = (*cmat)(r, c);
	mf->Close();

	//_corr = _corr.inverse();
	//std::cout << "_corr diag\n" << _corr.diagonal() << std::endl;

	if (kVerbosity) {
		std::cout << "Importing " << sysB-sysA
			  << " entries from " << mat_name
			  << " ( " << mat_file << " ) " << std::endl;
		std::cout << "Total number of systematics is "
			  << _corr.cols() << std::endl;
	}

	// define number of systematics!
	_nSys = _corr.cols();


	// now I am ready to import 1sigma histogram represeting linearisation of systematics
	_sysMatrix.clear();
	
	if (kVerbosity)
		std::cout << "BeamSample: creating matrix "
			  << _nBin << " x " << _nSys << std::endl;
	_sysMatrix[0] = Eigen::ArrayXXd::Zero(_nBin, _nSys);
	for (int sigma = -3; sigma < 4; sigma += 2)
		_sysMatrix[sigma] = Eigen::ArrayXXd::Zero(_nBin, _nSys);

	std::set<int> skip_sys;
	cd.Get("skip", skip_sys);	// errors to skip

	zeroEpsilons = true;
	//int off = 0;	// offset for global bin
	for (const std::string &it : _type) {
		// it is sample type name
		std::string file_name;
		if (!cd.Get("systematic_" + it, file_name)) {
			std::cout << "BeamSample: no systematic for " << it << " sample,"
				  << " parameters will be set to zero and not fitted\n";
			continue;
		}

		if (kVerbosity > 1)
			std::cout << "BeamSample: opening file " << file_name << std::endl;
		zeroEpsilons = false;

		TFile *sysf = new TFile(file_name.c_str());
		if (sysf->IsZombie()) {
			if (kVerbosity > 0)
				std::cerr << "WARNING - BeamSample: file " << file_name
					  << " does not exist" << std::endl;
			continue;
		}

		TIter next(sysf->GetListOfKeys());
		TH1D* hsys;
		TKey *k;
		while ((k = static_cast<TKey*>(next())))
		{
			std::string sysname = k->GetName();

			int sigma = 0, k_err = 0;
			if (sysname.find_first_of('_') != sysname.find_last_of('_')) {
			// spline systematics
				if (sysname.find("m3") != std::string::npos)
					sigma = -3;
				else if (sysname.find("m1") != std::string::npos)
					sigma = -1;
				else if (sysname.find("p1") != std::string::npos)
					sigma = 1;
				else if (sysname.find("p3") != std::string::npos)
					sigma = 3;

				hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));

				sysname.erase(sysname.find_last_of('_'));
				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number now
			}
			else {
				hsys = static_cast<TH1D*>(sysf->Get(sysname.c_str()));
				sysname.erase(0, sysname.find_first_of('_')+1);
				//sysname is just number
			}

			k_err = std::stoi(sysname);
			if (k_err < sysA || k_err >= sysB
			    || skip_sys.find(k_err) != skip_sys.end()) {
				if (kVerbosity)
					std::cout << "BeamSample: skipping " << hsys->GetTitle()
						<< " systematic" << std::endl;
				continue;
			}

			// past this point there should be something in the systematic file
			zeroEpsilons = false;

			int i = _offset[it];
			for (int n : _binpos[it]) {
				if (std::abs(hsys->GetBinContent(n+1) - 1) > 1)
					std::cout << k_err << ", " << n << " out of scale\n";
				if (sigma)	// it is a spline file, i.e. sigma != 0
					_sysMatrix[sigma](i, k_err - sysA)
						= sigma * (hsys->GetBinContent(n+1) - 1);
						//= sigma * std::min(1., std::max(-1., 
							   //hsys->GetBinContent(n+1) - 1));
				else		// not a spline, fill manually
					for (int s = -3; s < 4; s += 2)
						_sysMatrix[s](i, k_err - sysA)
							= s * (hsys->GetBinContent(n+1) - 1);
						//= s * std::min(1., std::max(-1., 
							   //hsys->GetBinContent(n+1) - 1));
				++i;
			}
		}
		sysf->Close();
	}
}

std::vector<std::pair<int, int> > BeamSample::AllSlices(std::string it, double skerr)
{
	if (!_nScale) {
		return Sample::AllSlices(it, skerr);
	}

	std::vector<std::pair<int, int> > allslices;

	// energy shift
	double shift = 1 + skerr * _scale[it];
	// offset between global bin and energy bin
	//int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// returns absolute bin position 
	int i = _offset[it];
	for (int n : _binpos[it]) {
		int m0 = std::max(StartingBin(it, shift, n), _binpos[it].front()) - n + i;
		int m1 = std::min(EndingBin(it, shift, n), _binpos[it].back()+1) - n + i;
		allslices.push_back(std::make_pair(m0, m1));
		++i;
	}

	// loop over bins of this sample
	//for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {
	//	int m0 = StartingBin(it, shift, n);
	//	int m1 = EndingBin(it, shift, n);
	//	allslices.push_back(std::make_pair(m0 + off, m1 + off));
	//}

	return allslices;
}


// FactorFn needs ChiSquared object to call
std::vector<Eigen::ArrayXd> BeamSample::AllScale(FactorFn factor, std::string it, double skerr)
{
	if (!_nScale) {
		return Sample::AllScale(factor, it, skerr);
	}

	std::vector<Eigen::ArrayXd> allfacts;

	// scale error value
	double scale_err = _scale[it];
	// energy shift
	double shift = 1 + skerr * scale_err;
	// offset between global bin and energy bin
	//int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// alias to binning
	std::vector<double> &global = _global[it];

	//// loop over bins of this sample
	//for (int n = _binpos[it].first; n < _binpos[it].second; ++n) {

		// scale systematic is related to bin number
		// scaling bin start
		//int m0 = StartingBin(it, shift, n);
		//int m1 = EndingBin(it, shift, n);

	for (int n : _binpos[it]) {
		// this is just relative bin position
		int m0 = std::max(StartingBin(it, shift, n), _binpos[it].front());
		int m1 = std::min(EndingBin(it, shift, n), _binpos[it].back()+1);
		//std::cout << "nonzero bin " << n << " is " << global[n]
		//	  << " - " << global[n+1]
		//	  << ", between " << m0 << " and " << m1 << std::endl;

		// unscaled/original bins
		//double b0_n = global[n - off];
		//double b1_n = global[n - off + 1];
		//
		// bin edges
		double b0_n = global[n];
		double b1_n = global[n + 1];

		std::vector<double> thisfact;
		for (int m = m0; m < m1; ++m) {

			double b0_m = global[m];
			double b1_m = std::min(global[m + 1], 30 / shift);
			//std::cout << "\tshift " << shift * b0_m << " - " << shift * b1_m << std::endl;
			// scaled bins
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
double Scale(FactorFn factor,
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
*/
