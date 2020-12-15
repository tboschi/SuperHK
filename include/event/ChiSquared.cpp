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

ChiSquared::ChiSquared(const std::string &card) :
	_nBin(-1),
	_nSys(-1)
{
	CardDealer cd(card);
	Init(cd);
}

ChiSquared::ChiSquared(const CardDealer &cd) :
	_nBin(-1),
	_nSys(-1)
{
	Init(cd);
}


ChiSquared::ChiSquared(CardDealer *cd) :
	_nBin(-1),
	_nSys(-1)
{
	Init(*cd);
}


void ChiSquared::Init(const CardDealer &cd)
{
	if (!cd.Get("lm_0", lm_0))
		lm_0 = 1;		//default value
	if (!cd.Get("lm_up", lm_up))
		lm_up = 5;		//default value
	if (!cd.Get("lm_down", lm_down))
		lm_down = 10;		//default value
	if (!cd.Get("lm_min", lm_min))
		lm_min = 0;		//default value

	if (!cd.Get("max_iterations", maxIteration))
		maxIteration = 10;
	if (!cd.Get("max_random_trials", maxTrials))
		maxTrials = 1e4;
	if (!cd.Get("fit_error", fitErr))
		fitErr = 1e-9;


	if (!cd.Get("verbose", kVerbosity))
		kVerbosity = 0;
}

bool ChiSquared::Combine()
{
	CombineBinning();
	CombineCorrelation();
	CombineSystematics();

	return (_nBin >= 0 && _nSys >= 0);
}


// total number of systematics
int ChiSquared::NumSys()
{
	return _nSys;
}


/*
int ChiSquared::NumBin() {
	// compute the number of nonempty bins
	// if not defined, the routine will do this check
	// and creates lower and upper bin pairs in _limits
	//
	if (_nBin < 0)
		_nBin = std::accumulate(_sample.begin(), _sample.end(), 0,
			[&](double sum, Sample* is)
			{ return sum + is->NumBin(); });

	return _nBin;
}
*/


/*
int ChiSquared::NumScale()
{
	if (_nScale < 0)
		_nScale = std::accumulate(_sample.begin(), _sample.end(), 0,
			[&](double sum, Sample* is)
			{ return sum + is->NumScale(); });

	return _nScale;
}
*/


int ChiSquared::DOF() {
	return std::abs(_nBin - _nSys);
}


// defines nBin too
void ChiSquared::CombineBinning() {
	_nBin = std::accumulate(_sample.begin(), _sample.end(), 0,
		[](double sum, std::shared_ptr<Sample> is)
		{ return sum + is->_nBin; });
}

// defines _nSys too
void ChiSquared::CombineCorrelation()
{
	if (_nSys < 0)
		_nSys = std::accumulate(_sample.begin(), _sample.end(), 0,
			[](double sum, std::shared_ptr<Sample> is)
			{ return sum + is->_nSys; });

	_corr = Eigen::MatrixXd::Identity(_nSys, _nSys);
	_nSys = 0;

	for (const auto &is : _sample) {
		int ns = is->_nSys;
		_corr.block(_nSys, _nSys, ns, ns) = is->_corr;
		_nSys += ns;
	}

	_corr = _corr.inverse();
}

void ChiSquared::CombineSystematics()
{
	zeroEpsilons = std::all_of(_sample.begin(), _sample.end(),
			[](std::shared_ptr<Sample> is) { return is->zeroEpsilons; } );

	_sysMatrix[0] = Eigen::ArrayXXd::Zero(_nBin, _nSys);
	for (int sigma = -3; sigma < 4; sigma += 2) {
		_sysMatrix[sigma] = Eigen::ArrayXXd::Zero(_nBin, _nSys);
		int row_off = 0, col_off = 0;
		for (const auto &is : _sample) {
			_sysMatrix[sigma].block(row_off, col_off, is->_nBin, is->_nSys) =
				is->_sysMatrix[sigma];
			row_off += is->_nBin;
			col_off += is->_nSys;
		}
	}
}

void ChiSquared::SetPoint(int p) {
	for (const auto &is : _sample)
		is->SetPoint(p);
}


std::map<std::string, Eigen::VectorXd> ChiSquared::BuildSamples(std::shared_ptr<Oscillator> osc) {
	std::map<std::string, Eigen::VectorXd> samples;
	for (const auto &is : _sample) {
		std::map<std::string, Eigen::VectorXd> sample = is->BuildSamples(osc);
		samples.insert(sample.begin(), sample.end());
	}

	return samples;
}

Eigen::VectorXd ChiSquared::ConstructSamples(std::shared_ptr<Oscillator> osc) {
	Eigen::VectorXd vect(_nBin);
	int off = 0;
	for (const auto &is : _sample) {
		vect.segment(off, is->_nBin) = is->ConstructSamples(osc);
		off += is->_nBin;
	}

	return vect;
}

//On is the true spectrum, En is the observed spectrum
//return time taken for computation
Eigen::VectorXd ChiSquared::FitX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En)
{
	//initialize epsil with zeroes

	Eigen::VectorXd epsil = Eigen::VectorXd::Zero(_nSys);
	double x2 = X2(On, En, epsil);

	if (std::isnan(x2)) {
		std::cerr << "ChiSquared: ERROR - X2 is nan - dumping vectors on stdout\n";
		std::cout << "ChiSquared: On " << On.transpose() << "\n\n";
		std::cout << "ChiSquared: En " << En.transpose() << "\n\n";
		throw std::logic_error("ChiSquared: starting X2 is nan and it shouldn't be\n");
	}

	if (zeroEpsilons)
		return epsil;

	Eigen::VectorXd best_eps = epsil, prev_eps = epsil;
	double best_x2 = x2;

	int tries = 0;
	double step = 1e-3;
	//double stepSize = 1.0/maxIteration;

	while (tries < maxIteration) {
		unsigned int code = MinimumX2(On, En, epsil, x2);
		std::cout << "Try (" << tries << ") : ";

		if (code == 0) {	// success!
			std::cout << "success!\n";
			return epsil;
		}
		else if (code == 1) {
			std::cout << "no convergence, dx2 " << best_x2-x2;
			if (best_x2 > x2) {
				//step = (best_eps - epsil).norm() / 2.0;
				//step /= lm_down;
				//step = std::min(1.0, std::abs(best_x2 - x2));
				//prev_eps = best_eps;
				step = std::min(step, std::abs(best_x2 - x2));
				best_eps = epsil;
				best_x2 = x2;
				std::cout << "-> new best! X2 = " << best_x2;
			}
			else
				step /= 10.;	// reset step size
			std::cout << std::endl;
		}
		else if (code == 2) {
			std::cout << "bad point ";
			if (best_x2 > x2) {
				//step = (best_eps - epsil).norm() / 2.0;
				//step /= lm_down;
				//step = std::min(1.0, std::abs(best_x2 - x2));
				//prev_eps = best_eps;
				best_eps = epsil;
				best_x2 = x2;
				std::cout << "-> new best! X2 = " << best_x2;
			}
			step /= 10.;	// reset step size
			//step *= lm_up;
			//step = 1;
			std::cout << std::endl;
		}

		// find a better point
		size_t rands = 0;
		do {
			epsil.setRandom();
			epsil = best_eps + epsil * step;
			x2 = X2(On, En, epsil);	//new initial value
			++rands;
		} while (x2 > best_x2 && rands < maxTrials);

		//if (step < fitErr)
		//	step = 0.1;
		//for (int i = 0; i < 50; ++i) {
			//std::cout << "* " << (epsil-best_eps).norm() << "\tfrom best "
				//" X2 is\t" << X2(On, En, epsil) << std::endl;
		//}


		//epsil.normalize();

		//epsil = prev_eps * (1 - step) + step * best_eps; // + epsil.normalized() * stepSize;
		std::cout << "new step (" << step << ") is " << (epsil-best_eps).norm()
			  << " from best after " << rands << " trials\n";
		//epsil = best_eps + step * epsil / epsil.norm();
		//epsil = step * epsil;
		//x2 = X2(On, En, epsil);	//new initial value

		++tries;
	}

	if (kVerbosity > 1)
		std::cout << "\nNumber of attempts: " << tries << std::endl;

	return best_eps;
}


//uses Levenberg-Marquardt algorithm
//quites when x2 variation falls below sens and if lambda < 1
//which means it iconvergin
unsigned int ChiSquared::MinimumX2(const Eigen::VectorXd &On,
			  const Eigen::VectorXd &En, 
			  Eigen::VectorXd &epsil,
			  double &x2)
			     //double alpha)
{
	int c = 0;	//counter
	double lambda = lm_0;

	Eigen::VectorXd delta = Eigen::VectorXd::Ones(_nSys);

	//x2 is initial x2 value
	//x2 = X2(On, En, epsil);
	
	double diff = 1;
	//double minj = 1;
	//while (diff > err && minj > NumSys() * err && minj > 1e-6 &&
	//       lambda < 1e6 && c < maxIteration)
	//while (delta.norm() > NumSys() * err && c < maxIteration)

	if (kVerbosity > 2)
		std::cout << "Minimising fit from x2: " << x2 << std::endl;

	// default lm_min = 0, cause if lambda == 0 risk of infinite loop
	while (lambda > lm_min && std::abs(diff / DOF()) > fitErr
	      && delta.norm() / _nSys > fitErr) {
	//while (std::abs(diff) > fitErr
	      //&& delta.norm() > fitErr) {
		++c;	//counter

		// build hessian and gradient/jacobian together
		// to save computation time
		Eigen::VectorXd jac(_nSys);
		Eigen::MatrixXd hes(_nSys, _nSys);
		JacobianHessian(jac, hes, On, En, epsil);

		//add diagonal to hesT hes
		double maxd = hes.diagonal().maxCoeff();
		hes.diagonal() += Eigen::VectorXd::Constant(_nSys, maxd * lambda);
		delta = hes.ldlt().solve(jac);
		if (delta.norm() > 100) {
			return 2;	// change step
		}

		//cos for uphill moves
		//double cosb = delta.dot(newdt) / (newdt.norm() * delta.norm());

		Eigen::VectorXd nextp = epsil - delta;	//next step
		//check if this step is good
		double obs_x2 = ObsX2(On, En, nextp);
		double sys_x2 = SysX2(nextp);

		// requring both terms positive also picks up nan value
		if (!(obs_x2 >= 0 && sys_x2 >= 0)) // bad points
			return 2;	//X2 cannot be computed

		diff = x2 - (obs_x2 + sys_x2);

		if (diff > 0) {	//next x2 is better, update lambda and epsilons
			lambda /= lm_down;
			epsil = nextp;
			x2 = obs_x2 + sys_x2;
		}
		else if (lambda > 0 && delta.norm() > 0)	//next x2 is worse
			lambda *= lm_up;	//nothing changes but lambda
						//if step is nonnull

		//std::cout << "\ndelta\n" << delta << std::endl;
		//if (false)
		if (kVerbosity > 2) {
			//std::cout << "epsil " << delta.transpose() << std::endl;
			std::cout << c << " -> l " << lambda
				  << ",\tstep: " << delta.norm() 
				  << ", X2: " << ObsX2(On, En, epsil)
				  << " + " << SysX2(epsil)
				  << " ( " << diff
				  << " ) " << diff / delta.norm() << std::endl;
		}
	}

	if (lambda <= lm_min)
		return 1;	// no convergence
	return (lambda <= lm_0) ? 0 : 1;
}


Eigen::VectorXd ChiSquared::Jacobian(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En, 
				      const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(_nSys);
	Eigen::MatrixXd hes(_nSys, _nSys);
	JacobianHessian(jac, hes, On, En, epsil);

	return 2 * jac;
}


Eigen::MatrixXd ChiSquared::Hessian(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En, 
				     const Eigen::VectorXd &epsil)
{
	Eigen::VectorXd jac(_nSys);
	Eigen::MatrixXd hes(_nSys, _nSys);
	JacobianHessian(jac, hes, On, En, epsil);

	return 2 * hes;
}


void ChiSquared::JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				 const Eigen::VectorXd &On,
				 const Eigen::VectorXd &En, 
				 const Eigen::VectorXd &epsil)
{
	//corr is inverse of correlation matrix
	jac = _corr * epsil;	//gradient/jacobian
	hes = _corr;		//hessian
	//hes.setZero();		//hessian

		//std::cout << "before jac\n" << jac.transpose() << "\n"
		//	  << "before hes\n" << hes.transpose() << std::endl;

	// event distribution with systematics applied
	const Eigen::ArrayXd Ep = Gamma(En, epsil);
	// matrix contains derivative terms for each systematics and for each bin
	const Eigen::ArrayXXd Fp = one_Fp(epsil);

	if (kVerbosity > 7)
		std::cout << "Obs : " << On.transpose() << "\n"
			  << "Exp : " << Ep.transpose() << std::endl;

	// loop over samples
	int err_off = 0, bin_off = 0;
	for (const auto &is : _sample) {	// beam or atmo
		for (const std::string &it: is->_type) {	// 1Re, 1Rmu, ...
			int t = is->ScaleError(it) + err_off;
			double skerr = is->_nScale ? epsil(t) : 0;

			std::vector<std::pair<int, int> > slices = is->AllSlices(it, skerr);
			std::vector<Eigen::ArrayXd> scales = is->AllScale(&Sample::Nor, it, skerr);
			std::vector<Eigen::ArrayXd> jacobs = is->AllScale(&Sample::Jac, it, skerr);
			std::vector<Eigen::ArrayXd> hessis = is->AllScale(&Sample::Hes, it, skerr);


	// looping over bins of given sample and type <is, it>
	/* >>>>>>>>>>>>>>> */
	for (size_t m = 0, n = is->_offset[it] + bin_off; m < is->_binpos[it].size(); ++m, ++n) {

		int m0 = slices[m].first, dm = slices[m].second - m0;
		m0 += bin_off;
		const Eigen::ArrayXd &ev = Ep.segment(m0, dm);

		if (kVerbosity > 4)
			std::cout << "type " << it << ", bin "
				  << is->_binpos[it][m] << " (abs " << n << ")"
				  << " scales are between "
				  << m0 << " and " << m0 + dm << std::endl;

		double on = On(n);
		double en = (scales[m] * ev).sum();

		double one_oe = 1 - on / en;

		double en_jac = (jacobs[m] * ev).sum();

		// jacobian
		if (is->_nScale) {
			jac(t) += one_oe * en_jac;
			if (kVerbosity > 5)
				std::cout << t << " E " << epsil(t) << ", jac "
					  << one_oe * en_jac << " ("
					  << jac(t) << ") = " << one_oe << ", "
					  << on << ", " << en << ", [ "
					  << scales[m].transpose() << " : "
					  << ev.transpose() << " ] "
					  << en_jac << std::endl;

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
				+ one_oe * (hessis[m] * ev).sum();
		}

		//do the rest, including mixed term with SK syst
		for (int k = err_off; k < err_off + is->_nSys - is->_nScale; ++k) {
			// scaled jacobian in this bin
			const Eigen::ArrayXd &kk = Fp.col(k).segment(m0, dm);

			// ev * kk is the derivative wrt to k-th syst
			double gn_jac = (scales[m] * ev * kk).sum();

			// jacobian
			jac(k) += one_oe * gn_jac;
			if (kVerbosity > 5) { // && kk.abs().sum() > 0) {
				//std::cout << "kk\t" << Fp.col(k).transpose() << std::endl;
				std::cout << k << " E " << epsil(k) << ", jac " << one_oe * gn_jac << " ("
					<< jac(k) << ") = " << one_oe << ", "
					<< on << ", " << en << ", [ "
					<< scales[m].transpose() << " : "
					<< ev.transpose() << " : "
					<< kk.transpose() << " ] " << gn_jac << std::endl;
			}

			//diagonal term
			hes(k, k) += on * pow(gn_jac / en, 2);

			if (is->_nScale)
				//mixed term with SK error
				hes(k, t) += en_jac * gn_jac * on / en / en
					+ one_oe * (jacobs[m] * ev * kk).sum();
					//+ one_oe * (jacobs[m] * kk).sum();

			// mixed terms with non energy scale parameters
			for (int j = k + 1; j < err_off + is->_nSys - is->_nScale; ++j) {

				if (kVerbosity > 6)
					std::cout << "m " << m << "\tk " << k
						  << "\tj " << j << "\n";
				const Eigen::ArrayXd &jj
					= Fp.col(j).segment(m0, dm);

				hes(k, j) += gn_jac * on / en / en
					  * (scales[m] * ev * jj).sum()
					+ one_oe
					* (scales[m] * ev * jj * kk).sum();
			}
		}
	}
	/* <<<<<<<<<<<<<<<<<<<<<< */

		}
		err_off += is->_nSys;
		bin_off += is->_nBin;
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

	return hes.inverse();
}

Eigen::VectorXd ChiSquared::Variance(const Eigen::VectorXd &On,
				     const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil)
{
	return Covariance(On, En, epsil).diagonal();
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

	int err_off = 0, bin_off = 0;
	for (const auto &is : _sample) {	// beam or atmo
		for (const std::string &it: is->_type) {	// 1Re, 1Rmu, ...
			if (kVerbosity > 4)
				std::cout << "type " << it << std::endl;
			int t = is->ScaleError(it) + err_off;
			double skerr = is->_nScale ? epsil(t) : 0;

			std::vector<std::pair<int, int> > slices = is->AllSlices(it, skerr);
			std::vector<Eigen::ArrayXd> scales = is->AllScale(&Sample::Nor, it, skerr);

			//for (int m = 0, n = is->_offset[it]; m < is->_binpos[it].size(); ++m, ++n) {
			for (size_t m = 0, n = is->_offset[it] + bin_off;
					   m < is->_binpos[it].size(); ++m, ++n) {
				//if (On(n) == 0)
				//	continue;
				int m0 = slices[m].first, dm = slices[m].second - m0;
				m0 += bin_off;
				double en = (scales[m] * gam.segment(m0, dm)).sum();
				chi2(n) = 2 * en - 2 * On(n) * (1 + log(en / On(n)));

				if (kVerbosity > 5)
					std::cout << "x2 @ " << m << ", " << n << " = " << chi2(n)
						  << " (" << en << ", " << En(n) << ", " << On(n) << ")" << std::endl;
			}
		}
		err_off += is->_nSys;
		bin_off += is->_nBin;
	}

	return chi2;
}

double ChiSquared::RawX2(const Eigen::VectorXd &On,
				 const Eigen::VectorXd &En)
{
	// modified expected events with systematics
	return RawX2n(On, En).sum();
}


Eigen::ArrayXd ChiSquared::RawX2n(const Eigen::VectorXd &On,
				  const Eigen::VectorXd &En)
{
	Eigen::ArrayXd chi2 = Eigen::ArrayXd::Zero(std::min(On.size(), En.size()));

	for (int n = 0; n < chi2.size(); ++n)
		if (On(n) > 0)
			chi2(n) = 2 * En(n) - 2 * On(n) * (1 + log(En(n) / On(n)));
	return chi2;
}

// return systematic X2
double ChiSquared::SysX2(const Eigen::VectorXd &epsil) {
	return epsil.transpose() * _corr * epsil;
}


// this return the spectrum modified with the systematics
Eigen::ArrayXd ChiSquared::Gamma(const Eigen::VectorXd &En,
				 const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXd gam = En.array();

	int err_off = 0;
	for (const auto &is : _sample) {	// beam or atmo
		for (int k = err_off; k < err_off + is->_nSys - is->_nScale; ++k) {
			auto oneFk = one_Fk(epsil(k), k);
			gam *= oneFk;//(epsil(k + off), k + off);
			//for (int n = 0; n < gam.size(); ++n)
			//	if (oneFk(n) < 0)
			//		std::cout << gam(n) << " with err " << k << " is " << epsil(k) << " at " << n << " : " << oneFk(n) << std::endl;

			//std::cout << k << ": " << one_Fk(epsil(k+off), k+off).transpose() << std::endl;
		}
		err_off += is->_nSys;
	}

	return gam;
}



// this is (1 + F)
Eigen::ArrayXXd ChiSquared::one_F(const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXXd f(_nBin, _nSys);
	
	int err_off = 0;
	for (const auto &is : _sample) {	// beam or atmo
		for (int k = err_off; k < err_off + is->_nSys - is->_nScale; ++k)
			f.col(k) = one_Fk(epsil(k), k);
		err_off += is->_nSys;
	}

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
	else {
		dl = 1; du = 3;
	}

	const Eigen::ArrayXd &sl = _sysMatrix[dl].col(k);
	Eigen::ArrayXd su = _sysMatrix[du].col(k);

	su = (su - sl) / (du - dl);
	return Eigen::ArrayXd::Ones(_nBin) + su * (err - dl) + sl;
}

// this is Fp / (1 + F) which is derivative wrt epsilon
Eigen::ArrayXXd ChiSquared::one_Fp(const Eigen::VectorXd &epsil)
{
	Eigen::ArrayXXd fp(_nBin, _nSys);
	
	int off = 0;
	for (const auto &is : _sample) {	// beam or atmo
		for (int k = 0; k < is->_nSys - is->_nScale; ++k)
			fp.col(k) = one_Fpk(epsil(k + off), k + off);
		off += is->_nSys;
	}

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
	else {
		dl = 1; du = 3;
	}

	Eigen::ArrayXd sl = _sysMatrix[dl].col(k);
	Eigen::ArrayXd su = _sysMatrix[du].col(k);

	su = (su - sl) / (du - dl);
	sl = Eigen::ArrayXd::Ones(_nBin) + su * (err - dl) + sl;

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
