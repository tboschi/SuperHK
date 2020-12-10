/*
 * Sample base class
 * Author: Tommaso Boschi
 */

#ifndef Sample_H
#define Sample_H

#include <iostream>
#include <string>
#include <map>
#include <memory>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

#include "Eigen/Dense"

class Sample
{
	public:
		Sample(const std::string &card) :
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

		Sample(const CardDealer &cd) :
			_nBin(-1),
			_nSys(-1),
			_nScale(0)
       		{
			if (!cd.Get("verbose", kVerbosity))
				kVerbosity = 0;
			if (!cd.Get("stats", _stats))
				_stats = 1.0;
		}

		Sample(CardDealer *cd) :
			_nBin(-1),
			_nSys(-1),
			_nScale(0)
       		{
			if (!cd->Get("verbose", kVerbosity))
				kVerbosity = 0;
			if (!cd->Get("stats", _stats))
				_stats = 1.0;
		}

		virtual ~Sample() = default;

		/*
		virtual int NumBin() {
			if (_nBin < 0) {
				std::throw("Sample : number of bins not defined");
			return _nBin;
		}

		virtual int NumSys() {
			if (_nSys < 0) {
				std::throw("Sample : number of systematics not defined");
			return _nSys;
		}

		virtual int NumScale() {
			if (_nScale < 0) {
				std::throw("Sample : number of scale systematics not defined");
			return _nScale;
		}
		*/


		virtual int ScaleError(std::string it) {
			return _nScale ? _nSys - _nScale + _type_scale[it] : 0;
		}

		// same for every one
		virtual void DefineBinning() {
			_nBin = 0;
			_allBin = 0;
			_offset.clear(); _binpos.clear();

			// unoscillated spectrum, just for counting nonempty bins
			std::map<std::string, Eigen::VectorXd> samples = BuildSamples();

			if (kVerbosity > 1)
				std::cout << "Sample: number of samples " << samples.size() << std::endl;
			for (const auto &is : samples) {
				// find first and last non empty bins
				//int b0 = 0, b1 = is.second.size();
				//for (int i = 1; i < is.second.size(); ++i)
				//	if (is.second(i) > 0) {	// first bin above zero
				//		b0 = i;
				//		break;
				//	}
				//for (int i = is.second.size()-1; i > 0; --i)
				//	if (is.second(i-1) > 0)	{ // last bin above zero
				//		b1 = i;
				//		break;
				//	}

				std::vector<size_t> bpos;
				bpos.reserve(is.second.size());

				// bpos has nonzero bins
				for (int i = 0; i < is.second.size(); ++i)
					if (is.second(i) > 0)
						bpos.push_back(size_t(i));

				// lims has start and end of nonzero bins
				//lims.push_back(_nBin + lims.front());
				//for (int b = 1; b < lims.size()-1; ++b)
				//	if (lims[b] - lims[b-1] > 1) {
				//		lims.push_back(_nBin + lims[b-1]);
				//		lims.push_back(_nBin + lims[b]);
				//	}
				//lims.push_back(_nBin + lims.back());

				if (kVerbosity > 2) {
					std::cout << "Sample: " << is.first << " has "
						  << bpos.size() << " nonzero bins";
					if (bpos.size())
						  std::cout << ", from " << bpos.front()
							  << " to " << bpos.back();
					std::cout << std::endl;
				}

				// store first and last non empty bins
				//_limits[is.first] = lims; //std::pair<int, int>(b0, b1);
				_offset[is.first] = _nBin; //std::pair<int, int>(b0, b1);

				// refers to position as global binning
				// eventually the samples will be flattened in one dimension,
				// so it is useful to know where the samples start and end
				_binpos[is.first] = bpos; //std::pair<int, int>(_nBin, _nBin + b1 - b0);

				// store global binning information for energy shift
				//_global[is.first] = is.second.rows();

				// count number of total bins
				_nBin += bpos.size();
				_allBin += is.second.size();

				//delete is.second;
			}

			if (kVerbosity)
				std::cout << "Sample: Number of nonzero bins: "
					<< _nBin << " / " << _allBin << std::endl;
		}

		// same for every one
		// BuildSpectrum and then collates everything on a Eigen::Vector
		virtual Eigen::VectorXd ConstructSamples(std::shared_ptr<Oscillator> osc = nullptr) {

			std::map<std::string, Eigen::VectorXd> samples = BuildSamples(osc);


			Eigen::VectorXd vect(_nBin);

			// compress all samples into one vector without zeros
			for (const auto &is : samples) {
				int i = _offset[is.first];
				for (int n : _binpos[is.first]) {
					vect(i) = is.second(n);
					++i;
				}
			}
				//const auto &lims = _limits[is.first];

				//int span = 0;
				//int prev = _limits[is.first].front();
				//for (int i : _limits[is.first]) {
				//	vect.segment(off, span) = is.second.segment(lims[i], span) * stats;
				//	off += span;
				//}
			//}

			return _stats * vect;
		}

		// must be defined in derived
		virtual void LoadReconstruction(const CardDealer &cd) = 0;
		virtual void LoadSystematics(const CardDealer &cd) = 0;

		// must be defined in derived
		virtual std::map<std::string, Eigen::VectorXd>
			BuildSamples(std::shared_ptr<Oscillator> osc = nullptr) = 0;


		// energy bin scaling routines

		/*
		virtual std::string TypeFromBin(int n) {
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
		*/

		virtual size_t StartingBin(std::string it, double shift, int n)
		{
			auto im = std::lower_bound(_global[it].begin(),
					_global[it].end(),
					_global[it][n] / shift);
			return std::distance(_global[it].begin(), im) - 1; // negative value?
			//return std::max(int(std::distance(_global[it].begin(), im)) - 1,
			//		_limits[it].first);
		}

		virtual size_t EndingBin(std::string it, double shift, int n)
		{
			auto im = std::upper_bound(_global[it].begin(),
					_global[it].end(),
					_global[it][n+1] / shift);
			return std::distance(_global[it].begin(), im);
			//return std::min(int(std::distance(_global[it].begin(), im)),
			//		_limits[it].second);
		}


		//function to compute the factor
		typedef double (Sample::*FactorFn)(const std::vector<double> &);


		// Apply energy scale dilation
		// FactorFn is a function which changes if calculating energy scale, jacobian or hessian
		// The dilation does not depend on En, therefore 
		// En can be anything

		virtual std::vector<std::pair<int, int> > AllSlices(std::string it, double skerr) {

			std::vector<std::pair<int, int> > allslices;

			// energy shift
			//double shift = 1 + skerr * _scale[it];
			// offset between global bin and energy bin
			//int off = _binpos[it].first - _limits[it].first;
			// absolute systematic error for this scale parameter

			// loop over bins of this sample
			//for (int n = _binpos[it].first; n < _binpos[it].second; ++n)
			for (size_t n = _offset[it]; n < _offset[it] + _binpos[it].size(); ++n)
				allslices.push_back(std::make_pair(n, n+1));

			return allslices;
		}


		virtual std::vector<Eigen::ArrayXd> AllScale(FactorFn factor,
				std::string it, double skerr)
		{
			std::vector<Eigen::ArrayXd> allfacts;
			//for (int n = _binpos[it].first; n < _binpos[it].second; ++n)
			for (size_t n = _offset[it]; n < _offset[it] + _binpos[it].size(); ++n)
				allfacts.push_back(Eigen::ArrayXd::Ones(1));
			return allfacts;
		}


		// derivative terms for chi2
		virtual double Nor(const std::vector<double> &term) {
			//return f / shift / db;
			return term[2] / term[1] / term[4];
		}

		virtual double Jac(const std::vector<double> &term) {
			//return scale_err / shift / db * (fd - f / shift)
			return term[0] / term[1] / term[4] * (term[3] - term[2] / term[1]);
		}

		virtual double Hes(const std::vector<double> &term) {
			//return pow(scale_err / shift, 2) / db * (f / shift - fd)
			return 2 * pow(term[0] / term[1], 2) / term[4]
				//* (term[2] / term[1] - term[3]);
				* (term[2] / term[1] - 2 * term[3]);
		}


		void SetPoint(int p) {
			_point = p;
		}


	protected:
		int kVerbosity;
		bool zeroEpsilons;

		// matter profile stored here
		// remember to set it with
		// 	osc->SetMatterProfile(_lens_dens)
		Oscillator::Profile _lens_dens;

		std::set<std::string> _type;
		std::map<std::string, std::pair<Nu::Flavour, Nu::Flavour> > _oscf;

		// for binning information
		int _nBin, _allBin;
		//std::map<std::string, std::pair<int, int> > _binpos;
		//std::map<std::string, std::pair<int, int> > _limits;
		std::map<std::string, std::vector<size_t> > _binpos;
		std::map<std::string, size_t> _offset;
		std::map<std::string, std::vector<double> > _global;
		// store point for pre computed bins
		int _point;
		double _stats;

		// for systematics
		int _nSys;
		Eigen::MatrixXd _corr;
		std::map<int, Eigen::ArrayXXd> _sysMatrix;

		// for scale errors
		int _nScale;
		std::map<std::string, double> _scale;
		std::map<std::string, int> _type_scale;

		// ChiSquared is friend and so can access private members
		friend class ChiSquared;
};

#endif

/*
		virtual double Scale(FactorFn factor,
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
