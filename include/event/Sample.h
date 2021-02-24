/*
 * Sample base class
 * Author: Tommaso Boschi
 */

#ifndef Sample_H
#define Sample_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <memory>
#include <utility>

#include "tools/CardDealer.h"
#include "physics/Const.h"
#include "physics/Flavors.h"
#include "physics/Oscillator.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TKey.h"
#include "TMatrixD.h"
#include "TMatrixT.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

class Sample
{
	public:
		Sample(const std::string &card);
		Sample(const CardDealer &cd);
		Sample(CardDealer *cd);

		virtual ~Sample() = default;

		virtual void Load(const CardDealer &cd, std::string process = "");

		// must be defined in derived class!!
		virtual void LoadReconstruction(const CardDealer &cd) = 0;
		virtual void LoadSystematics(const CardDealer &cd) = 0;
		virtual std::unordered_map<std::string, Eigen::VectorXd>
			BuildSamples(std::shared_ptr<Oscillator> osc = nullptr) = 0;

		virtual int NumBin();
		virtual int NumSys();
		virtual int NumScale();
		virtual std::unordered_map<std::string, size_t> ScaleError();
		virtual int ScaleError(std::string it);

		// same for every one
		virtual void DefineBinning(bool zerosuppress = true);

		// same for every one
		// BuildSpectrum and then collates everything on a Eigen::Vector
		virtual Eigen::VectorXd ConstructSamples(std::shared_ptr<Oscillator> osc = nullptr);
		// opposite function as above but must be defined in child class
		virtual std::unordered_map<std::string, Eigen::VectorXd>
			Unfold(const Eigen::VectorXd &En) = 0;


		// return bins, must be impemented in base class
		virtual std::vector<double> GetErecoBins(std::string it);
		virtual std::vector<double> GetEtrueBins(std::string it);

		virtual double X2(const Eigen::VectorXd &On,
			  const Eigen::VectorXd &En,
			  const Eigen::VectorXd &epsil);
		virtual double SysX2(const Eigen::VectorXd &epsil);
		virtual double ObsX2(const Eigen::VectorXd &On,
			     const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);
		virtual Eigen::ArrayXd ObsX2n(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);

		// method available to anybody to compute just a raw chi2, no systematics
		static double RawX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En) {
			return RawX2n(On, En).sum();
		}

		static Eigen::ArrayXd RawX2n(const Eigen::VectorXd &On, const Eigen::VectorXd &En) {
			assert((On.size() == En.size()) && "ChiSquared: RawX2n On and En different sizes");

			const Eigen::ArrayXd &on = On.array();
			const Eigen::ArrayXd &en = En.array();

			Eigen::ArrayXd chi2 = 2 * en - 2 * on * (1 + en.log() - on.log());
			//return chi2;
			return (chi2.isFinite()).select(chi2, 0);
		}


		// apply systematics of this sample
		virtual Eigen::VectorXd Gamma(const Eigen::VectorXd &En, const Eigen::VectorXd &epsil);
		virtual Eigen::VectorXd GammaP(const Eigen::VectorXd &En, const Eigen::VectorXd &epsil);

		virtual Eigen::ArrayXd one_Fk(double err, int k);
		virtual Eigen::ArrayXd one_Fpk(double err, int k);

		virtual Eigen::ArrayXXd one_F(const Eigen::VectorXd &epsil);
		virtual Eigen::ArrayXXd one_Fp(const Eigen::VectorXd &epsil);

		// energy scale stuff
		//function to compute the xi factor (please refer to the doc)
		using Xi = double (Sample::*)(const std::array<double, 5> &);

		//virtual size_t StartingBin(std::string it, double shift, int n);
		//virtual size_t EndingBin(std::string it, double shift, int n);
		//virtual std::vector<std::pair<int, int> > AllSlices(std::string it, double skerr);
		//virtual std::vector<Eigen::ArrayXd> AllScale(Xi xi, std::string it, double skerr);
		virtual double Nor(const std::array<double, 5> &term);
		virtual double Jac(const std::array<double, 5> &term);
		virtual double Hes(const std::array<double, 5> &term);

		virtual Eigen::SparseMatrix<double> ScaleMatrix(Xi xi, const Eigen::VectorXd &epsil);

	protected:
		int kVerbosity;
		bool zeroEpsilons;

		// matter profile stored here
		// remember to set it with
		// 	osc->SetMatterProfile(_lens_dens)
		Oscillator::Profile _lens_dens;

		std::set<std::string> _type;
		std::unordered_map<std::string, std::pair<Nu::Flavor, Nu::Flavor> > _oscf;

		// for binning information
		int _nBin; //, _allBin;
		std::unordered_map<std::string, size_t> _offset;
		std::unordered_map<std::string, std::vector<size_t> > _binpos;
		std::unordered_map<std::string, std::vector<double> > _global_true, _global_reco;
		// store point for pre computed bins
		size_t _point;
		double _stats;

		// for systematics
		int _nSys;
		Eigen::MatrixXd _corr;
		std::unordered_map<int, Eigen::ArrayXXd> _sysMatrix;

		// for scale errors
		int _nScale;
		std::unordered_map<std::string, std::pair<double, size_t> > _scale;

		// ChiSquared is friend and so can access private members
		friend class ChiSquared;
};

#endif

/*
 * old Scaling function
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
