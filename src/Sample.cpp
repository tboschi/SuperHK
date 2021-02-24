#include "event/Sample.h"


Sample::Sample(const std::string &card) :
	_nBin(-1),
	_nSys(-1),
	_nScale(0)
{
	CardDealer cd(card);

	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;
	if (!cd.Get("stats", _stats))
		_stats = 1.0;
}

Sample::Sample(const CardDealer &cd) :
	_nBin(-1),
	_nSys(-1),
	_nScale(0)
{
	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;
	if (!cd.Get("stats", _stats))
		_stats = 1.0;
}

Sample::Sample(CardDealer *cd) :
	_nBin(-1),
	_nSys(-1),
	_nScale(0)
{
	if (!cd->Get("verbose", kVerbosity))
		kVerbosity = 0;
	if (!cd->Get("stats", _stats))
		_stats = 1.0;
}

void Sample::Load(const CardDealer &cd, std::string process) {
	// by default load reco, define bins without empty ones, load systs
	if (process.empty())
		process = "RBSE";
	else	// make everything upper
		std::transform(process.begin(), process.end(), process.begin(),
				[](unsigned char c){ return std::toupper(c); });

	if (process.find('R') != std::string::npos)
		LoadReconstruction(cd);
	if (process.find('B') != std::string::npos)
		// if no E option, keep also empty bins
		DefineBinning(process.find('E') != std::string::npos);
	if (process.find('S') != std::string::npos)
		LoadSystematics(cd);
}

int Sample::NumBin() {
	if (_nBin < 0)
		throw std::logic_error("Sample : number of bins not defined");
	return _nBin;
}

int Sample::NumSys() {
	if (_nSys < 0)
		throw std::logic_error("Sample : number of systematics not defined");
	return _nSys;
}

int Sample::NumScale() {
	if (_nScale < 0)
		throw std::logic_error("Sample : number of scale systematics not defined");
	return _nScale;
}

std::unordered_map<std::string, size_t> Sample::ScaleError() {
	std::unordered_map<std::string, size_t> scales;
	for (const auto & s : _scale)
		scales.emplace(s.first, s.second.second);
	return scales;
}

int Sample::ScaleError(std::string it) {
	if (_type.find(it) == _type.end())
		throw std::invalid_argument("Sample: unknown sample type " + it);
	return _nScale ? _scale[it].second : 0;
}

// common to all samples
void Sample::DefineBinning(bool zerosuppress) {
	_nBin = 0;
	_offset.clear(); _binpos.clear();

	// unoscillated spectrum, just for counting nonempty bins
	std::unordered_map<std::string, Eigen::VectorXd> samples = BuildSamples();

	if (kVerbosity > 1) {
		std::cout << "Sample: number of samples " << samples.size() << std::endl;
		if (zerosuppress)
			std::cout << "Sample: removing empty bins\n";
	}
	for (const std::string it : _type) {
		std::vector<size_t> bpos;
		bpos.reserve(samples[it].size());

		// bpos has only nonzero bins if zerosuppress is true
		for (int i = 0; i < samples[it].size(); ++i)
			if (!zerosuppress || samples[it](i) > 0)
				bpos.emplace_back(i);

		if (kVerbosity > 2) {
			std::cout << "Sample: " << it << " has "
				<< bpos.size() << " nonzero bins";
			if (bpos.size())
				std::cout << ", from " << bpos.front()
					<< " to " << bpos.back();
			std::cout << std::endl;
		}


		// refers to position as global binning
		// eventually the samples will be flattened in one dimension,
		// so it is useful to know where the samples start and end
	
		//_allBin[it] += samples[it].size();
		_offset[it] = _nBin;

		_nBin += bpos.size();
		_binpos[it] = std::move(bpos);
	}

	if (kVerbosity)
		std::cout << "Sample: Number of nonzero bins: " << _nBin << "\n";
}

// compress all samples into one vector without zeros
Eigen::VectorXd Sample::ConstructSamples(std::shared_ptr<Oscillator> osc) {

	// build samples return the full spectrum (zero bins included)
	std::unordered_map<std::string, Eigen::VectorXd> samples = BuildSamples(osc);
	Eigen::VectorXd vect(_nBin);

	for (const std::string it : _type) {
		int i = _offset[it];
		for (int n : _binpos[it]) {
			vect(i) = samples[it](n);
			++i;
		}
	}

	// rescale statistics -> move _stats inside buildsamples
	return _stats * vect;
}

std::vector<double> Sample::GetErecoBins(std::string it) {
	if (_type.find(it) == _type.end())
		throw std::invalid_argument("Sample: unknown sample type " + it);
	return _global_reco[it];
}

std::vector<double> Sample::GetEtrueBins(std::string it) {
	if (_type.find(it) == _type.end())
		throw std::invalid_argument("Sample: unknown sample type " + it);
	return _global_true[it];
}


double Sample::X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		  const Eigen::VectorXd &epsil)
{
	return ObsX2(On, En, epsil) + SysX2(epsil);
}


// epsil is an array with all sigmas including SK energy scale
// sum up all bin contributions to the chi2
double Sample::ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			 const Eigen::VectorXd &epsil)
{
	// modified expected events with systematics
	return ObsX2n(On, En, epsil).sum();
}


// return statistics X2 as a vector
// so X2 contribution to each non-zer bin
Eigen::ArrayXd Sample::ObsX2n(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			      const Eigen::VectorXd &epsil)
{
	assert((On.size() == _nBin));

	// modified expected events with systematics
	//return Sample::RawX2n(On, GammaP(En, epsil));
	auto raw = RawX2n(On, GammaP(En, epsil));
	//std::cout << "raw " << raw.sum() << " is " << raw.transpose() << "\n";
	return raw;
}

// return systematic X2
double Sample::SysX2(const Eigen::VectorXd &epsil) {
	return epsil.transpose() * _corr * epsil;
}

// this return the spectrum modified with the systematics
Eigen::VectorXd Sample::Gamma(const Eigen::VectorXd &En,
				 const Eigen::VectorXd &epsil)
{
	assert((En.size() == _nBin) && "Sample: not right number of entries");
	assert((epsil.size() == _nSys) && "Sample: not right number of systematic errors");

	Eigen::ArrayXd gam = En.array();		
	for (int k = 0; k < _nSys - _nScale; ++k)
		gam *= one_Fk(epsil(k), k);

	return gam.matrix();
}

// calling gamma function on the correct piece and apply scale systematic if any
// N.B. the scale systematics does not modify zero compression of spectrum
// and only nonzero bins are scaled, meaning any spillover is truncated
Eigen::VectorXd Sample::GammaP(const Eigen::VectorXd &En, const Eigen::VectorXd &epsil)
					 //std::shared_ptr<Oscillator> osc = nullptr)
{
	/*
	Eigen::VectorXd gam = Gamma(En, epsil);

	for (const std::string it : _type) {
		double skerr = _nScale ? epsil(_scale[it].second) : 0;

		std::vector<std::pair<int, int> > slices = AllSlices(it, skerr);
		std::vector<Eigen::ArrayXd> scales = AllScale(&Sample::Nor, it, skerr);
		for (size_t m = _offset[it]; m < _binpos[it].size(); ++m) {
			int m0 = slices[m].first, dm = slices[m].second - m0;
			if (dm == 0)
				gam(m) = 0;
			else
				gam(m) = (scales[m] * gam.segment(m0, dm).array()).sum();
		}
	}
	return gam;
	*/

	return ScaleMatrix(&Sample::Nor, epsil) * Gamma(En, epsil);
}


// Terms including systematic errors

// this is (1 + F) as matrix
Eigen::ArrayXXd Sample::one_F(const Eigen::VectorXd &epsil)
{
	// make sure right number of errors is passed
	assert((epsil.size() == _nSys) && "Sample: one_F epsilon passed wrong number of entries");

	Eigen::ArrayXXd f(_nBin, _nSys);
	
	for (int k = 0; k < _nSys - _nScale; ++k)
		f.col(k) = one_Fk(epsil(k), k);

	return f;
}

// k-th entry of the previous function
Eigen::ArrayXd Sample::one_Fk(double err, int k)
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
	else {
		dl = 1; du = 3;
	}

	const Eigen::ArrayXd &sl = _sysMatrix[dl].col(k);
	const Eigen::ArrayXd &su = _sysMatrix[du].col(k);

	//su = (su - sl) / (du - dl);
	//return (1 + su * (err - dl) + sl);
	return (1 + (su - sl) * (err - dl) / (du - dl) + sl);
}

// this is Fp / (1 + F) which is derivative wrt epsilon
Eigen::ArrayXXd Sample::one_Fp(const Eigen::VectorXd &epsil)
{
	// make sure right number of errors is passed
	assert((epsil.size() == _nSys) && "Sample: one_Fp epsilon passed wrong number of entries");

	Eigen::ArrayXXd fp(_nBin, _nSys);
	
	for (int k = 0; k < _nSys - _nScale; ++k)
		fp.col(k) = one_Fpk(epsil(k), k);

	return fp;
}

// k-th entry of the previous function
Eigen::ArrayXd Sample::one_Fpk(double err, int k)
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
	else {
		dl = 1; du = 3;
	}

	const Eigen::ArrayXd &sl = _sysMatrix[dl].col(k);
	const Eigen::ArrayXd &su = _sysMatrix[du].col(k);

	//su = (su - sl) / (du - dl);
	//sl = 1 + su * (err - dl) + sl;
	//return su / sl;

	return (su - sl) / ((1 + sl) * (du - dl) + (su - sl) * (err - dl));
}


Eigen::SparseMatrix<double> Sample::ScaleMatrix(Xi xi, const Eigen::VectorXd &epsil)
{
	Eigen::SparseMatrix<double> scale(_nBin, _nBin);
	scale.reserve(Eigen::VectorXi::Constant(_nBin, 1));
	for (int i = 0; i < _nBin; ++i)
		scale.insert(i, i) = 1;
	//scale.makeCompressed();
	return scale;
}


// derivative terms for chi2
double Sample::Nor(const std::array<double, 5> &term) {
	//return f / shift / db;
	return term[2] / term[1] / term[4];
}

double Sample::Jac(const std::array<double, 5> &term) {
	//return scale_err / shift / db * (fd - f / shift)
	return term[0] / term[1] / term[4] * (term[3] - term[2] / term[1]);
}

double Sample::Hes(const std::array<double, 5> &term) {
	//return pow(scale_err / shift, 2) / db * (f / shift - fd)
	return 2 * pow(term[0] / term[1], 2) / term[4]
		//* (term[2] / term[1] - term[3]);
		* (term[2] / term[1] - 2 * term[3]);
}

/*
size_t Sample::StartingBin(std::string it, double shift, int n)
{
	auto im = std::lower_bound(_global_reco[it].begin(),
			_global_reco[it].end(),
			_global_reco[it][n] / shift);
	return std::max(int(std::distance(_global_reco[it].begin(), im)) - 1, 0); // negative value?
}

size_t Sample::EndingBin(std::string it, double shift, int n)
{
	auto im = std::upper_bound(_global_reco[it].begin(),
			_global_reco[it].end(),
			_global_reco[it][n+1] / shift);
	return std::distance(_global_reco[it].begin(), im);
}
*/


// Apply energy scale dilation
// is a function which changes if calculating energy scale, jacobian or hessian
// The dilation does not depend on En, therefore 
// En can be anything

/*
std::vector<std::pair<int, int> > Sample::AllSlices(std::string it, double skerr) {

	std::vector<std::pair<int, int> > allslices;

	// energy shift
	//double shift = 1 + skerr * _scale[it];
	// offset between local bin and energy bin
	//int off = _binpos[it].first - _limits[it].first;
	// absolute systematic error for this scale parameter

	// loop over bins of this sample
	for (size_t n = _offset[it]; n < _offset[it] + _binpos[it].size(); ++n)
		allslices.push_back(std::make_pair(n, n+1));

	return allslices;
}


std::vector<Eigen::ArrayXd> Sample::AllScale(Xi xi, std::string it, double skerr)
{
	std::vector<Eigen::ArrayXd> allfacts;
	for (size_t n = _offset[it]; n < _offset[it] + _binpos[it].size(); ++n)
		allfacts.push_back(Eigen::ArrayXd::Ones(1));
	return allfacts;
}
*/
