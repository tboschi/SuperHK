#include "event/BeamSample.h"

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
	//std::map<std::string, std::pair<flavor, flavor> >
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

	std::map<std::string, double> scale;
	if (!cd.Get("scale_error_", scale)) {
		double err;
		if (cd.Get("scale_error", err)) {
			for (const std::string &it : _type)
				_scale.emplace(it, std::make_pair(err, 0));
			_nScale = 1;
		}
		else
			_nScale = 0;
	}
	else { // length of scale is number of errors
		_nScale = 0;
		for (const auto &s : scale) {
			for (const std::string &it : _type) {
				if (it.find(s.first) != std::string::npos)
					_scale.emplace(it, std::make_pair(s.second, _nScale));
				else
					throw std::invalid_argument("BeamSample: not accepted string for scale error \"" + s.first + "\"");
			}
			++_nScale;	// <- contains right number
		}
	}



	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;

	if (kVerbosity) {
		std::cout << "BeamSample: types: ";
		for (const auto &it : _type)
			std::cout << "\t" << it;
		std::cout << "\nBeamSample: number of scale systematics " << _nScale << "\n";
		for (const auto &it : _scale)
			std::cout << "\t" << it.first << " -> " << it.second.first << "\n";
	}

	std::string profile;
	if (!cd.Get("density_profile", profile) && kVerbosity)
		std::cout << "BeamSample: density profile not set in card, please set it yourself" << "\n";
	else
		_lens_dens = Oscillator::GetMatterProfile(profile);

	// dispatch method from base class
	Load(cd, process);
}


void BeamSample::LoadReconstruction(const CardDealer &cd)
{
	// open it if already exists
	std::vector<std::string> reco_files;
	if (!cd.Get("reco_input", reco_files))
		throw std::invalid_argument("BeamSample: no reconstruction files in card,"
					    "very bad!");
	for (const std::string &reco : reco_files){
		LoadReconstruction(reco);
	}

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

	std::unordered_map<std::string, std::string> mSD;
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
			const double *by = h2->GetYaxis()->GetXbins()->GetArray();

			if (!_global_true.count(type))  {
				_global_true[type].assign(bx, bx + xs + 1);
				_global_reco[type].assign(by, by + ys + 1);
			}
			//_binX[is.first].assign(bx, bx + xs + 1);
			//_binY[is.first].assign(by, by + ys + 1);
		}

	inFile->Close();
}


std::unordered_map<std::string, Eigen::VectorXd>
	BeamSample::BuildSamples(std::shared_ptr<Oscillator> osc)
{
	if (osc)
		osc->SetMatterProfile(_lens_dens);

	std::unordered_map<std::string, Eigen::VectorXd> samples;
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
			const auto &bins = _global_true[hname]; //type?

			for (size_t i = 0; i < bins.size() - 1; ++i) {
				probs(i) = osc->Probability(chan.first, chan.second,
						(bins[i] + bins[i+1]) / 2.);
			}
			//std::cout<< "Number of bins " << bins.size() <<std::endl;
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
				if (kVerbosity && std::abs(hsys->GetBinContent(n+1) - 1) > 1)
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

	// update scale systematic position
	for (auto &s : _scale)
		s.second.second += _nSys - _nScale;
}

// decompress spectrum vector
std::unordered_map<std::string, Eigen::VectorXd> BeamSample::Unfold(const Eigen::VectorXd &En)
{
	assert(En.size() == _nBin);
	std::unordered_map<std::string, Eigen::VectorXd> samples;
	for (const std::string it : _type) {
		samples[it] = Eigen::VectorXd::Zero(_global_true[it].size() - 1);
		for (size_t m = 0, n = _offset[it]; m < _binpos[it].size(); ++m, ++n)
			samples[it](_binpos[it][m]) = En(n);
	}

	return samples;
}


Eigen::SparseMatrix<double> BeamSample::ScaleMatrix(Xi xi, const Eigen::VectorXd &epsil)
{
	if (!_nScale) {
		return Sample::ScaleMatrix(xi, epsil);
	}

	// n, m matrix 
	Eigen::SparseMatrix<double> scale(_nBin, _nBin);
	scale.reserve(Eigen::VectorXi::Constant(_nBin, 5));

	for (const std::string it : _type) {
		// scale error value
		double skerr = epsil(_scale[it].second);
		double shift = 1 + skerr * _scale[it].first;

		// alias to reco binning
		std::vector<double> &reco = _global_reco[it];

		int i = _offset[it];
		for (size_t n : _binpos[it]) {

			// this is just relative bin position
			//int m0 = std::max(StartingBin(it, shift, n), _binpos[it].front());
			//int m1 = std::min(EndingBin(it, shift, n), _binpos[it].back()+1);
			//
			// unscaled/original bins
			//double b0_n = local[n - off];
			//double b1_n = local[n - off + 1];
			//
			// bin edges
			double b0_n = reco[n];
			double b1_n = reco[n + 1];

			auto im = std::lower_bound(reco.begin(), reco.end(), b0_n / shift);
			size_t m0 = std::max(size_t(std::distance(reco.begin(), im)),
					  _binpos[it][0] + 1) - 1;

			for (size_t m = m0; m < reco.size() - 1; ++m) {
				double b0_m = reco[m];
				double b1_m = std::min(reco[m + 1], 30 / shift);
				// scaled bins
				//continue/break computation because there is no overlap
				if (shift * b0_m > b1_n)
					break;

				double f = (b1_n - b0_n + shift * (b1_m - b0_m)
						- std::abs(b0_n - shift * b0_m)
						- std::abs(b1_n - shift * b1_m)) / 2.;

				int s0 = b0_n - shift * b0_m < 0 ? -1 : 1;
				int s1 = b1_n - shift * b1_m < 0 ? -1 : 1;
				double ss = f > 0 ? 1 : 0.5;	//continuity factor
				double fd = ss * (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2.;

				std::array<double, 5> terms{{_scale[it].first, shift,
							     f, fd, b1_m - b0_m}};
				//std::cout << " at " << i << ", " << m - n + i << " = " << (this->*xi)(terms) << "\n";
				scale.insert(i, m - n + i) = (this->*xi)(terms);
			}
			++i;
		}
	}
	// not necessary to compress!
	//scale.makeCompressed();

	return scale;
}

/*
std::vector<std::pair<int, int> > BeamSample::AllSlices(std::string it, double skerr)
{
	if (!_nScale) {
		return Sample::AllSlices(it, skerr);
	}

	std::vector<std::pair<int, int> > allslices;

	// energy shift
	double shift = 1 + skerr * _scale[it].first;
	std::cout << "shift is " << shift << std::endl;
	std::cout << "skerr is " << skerr << std::endl;
	std::cout << "_scale[it] is " << _scale[it] << std::endl;
		
	// offset between global bin and energy bin
	//int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// returns absolute bin position
	//std::cout << "In all slices " <<std::endl; 
	int i = _offset[it];
	allslices.reserve(_binpos[it].size());
	for (int n : _binpos[it]) {
		std::cout << "n is " << n << std::endl;
		std::cout << "i is " << i << std::endl;
		std::cout << "StartingBin is " << StartingBin(it, shift, n) << std::endl;
		std::cout << "EndingBin is " << EndingBin(it, shift, n) << std::endl;
		std::cout << "_binpos[it].front() is " << _binpos[it].front() << std::endl;
		std::cout << "_binpos[it].back()+1 is " << _binpos[it].back()+1 << std::endl;
		
		int m0 = std::max(StartingBin(it, shift, n), _binpos[it].front());
		int m1 = std::min(EndingBin(it, shift, n), _binpos[it].back()+1);
		if (m0 > m1)
			m0 = m1;
		allslices.emplace_back(m0 - n + i, m1 - n + i);
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
std::vector<Eigen::ArrayXd> BeamSample::AllScale(Xi xi, std::string it, double skerr)
{
	if (!_nScale) {
		return Sample::AllScale(xi, it, skerr);
	}

	std::vector<Eigen::ArrayXd> allfacts;

	// scale error value
	double scale_err = _scale[it].first;
	// energy shift
	double shift = 1 + skerr * scale_err;
	// offset between global bin and energy bin
	//int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// alias to reco binning
	std::vector<double> &reco = _global_reco[it];

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
		//std::cout << "nonzero bin " << n << " is " << local[n]
		//	  << " - " << local[n+1]
		//	  << ", between " << m0 << " and " << m1 << std::endl;

		// unscaled/original bins
		//double b0_n = local[n - off];
		//double b1_n = local[n - off + 1];
		//
		// bin edges
		double b0_n = reco[n];
		double b1_n = reco[n + 1];

		std::vector<double> thisfact;
		for (int m = m0; m < m1; ++m) {

			double b0_m = reco[m];
			double b1_m = std::min(reco[m + 1], 30 / shift);
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

			std::array<double, 5> terms{{scale_err, shift, f, fd, b1_m - b0_m}};
			thisfact.push_back((this->*xi)(terms));
		}

		Eigen::ArrayXd arr = Eigen::Map<Eigen::ArrayXd>
			(thisfact.data(), thisfact.size(), 1);
		allfacts.push_back(arr);
	}

	return allfacts;
}
*/

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
	double b0_n = _local[it][n - off];
	double b1_n = _local[it][n - off + 1];

	//std::cout << "type " << it << " from " << _binpos[it].first << " to " << _binpos[it].second << std::endl;
	//std::cout << "bin " << n << " : " << b0_n << " - " << b1_n << std::endl;
	//std::cout << "starting from " << m0 << ": "
	//<< shift * _global[it][m0] << " and " << shift * _global[it][m0 + 1] << "\n";
	//std::cout << "off_global

	double ret = 0;
	// Loop over unscaled/original edges only nonempty bins
	for (int m = m0; m < _local[it].size()-1; ++m) {
		// scaled bins
		double b0_m = _local[it][m];
		double b1_m = std::min(_local[it][m + 1], 30 / shift);
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

