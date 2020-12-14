/* ChiSquared
 * compute the likelihood ratio chi-squared using the pull approach to include systematics
 */

#ifndef ChiSquared_H
#define ChiSquared_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <utility>
#include <numeric>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixT.h"
#include "TKey.h"

#include "Eigen/Dense"
#include "Eigen/SVD"

#include "tools/CardDealer.h"
#include "event/Sample.h"
#include "event/AtmoSample.h"
#include "event/BeamSample.h"
#include "physics/Oscillator.h"

class ChiSquared
{
	public:
		ChiSquared(const std::string &card);
		ChiSquared(const CardDealer &cd);
		ChiSquared(CardDealer *cd);

		template <class S>
		void Add(std::string card) {
			_sample.push_back(std::shared_ptr<Sample>(new S(card)));
		}
		void SetPoint(int p);

		void Init(const CardDealer &cd);
		bool Combine();
		void CombineBinning();
		void CombineCorrelation();
		void CombineSystematics();
		std::map<std::string, Eigen::VectorXd> BuildSamples(std::shared_ptr<Oscillator> osc = nullptr);
		Eigen::VectorXd ConstructSamples(std::shared_ptr<Oscillator> osc = nullptr);

		int NumSys();
		int NumBin();
		int DOF();

		Eigen::VectorXd FitX2(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En);
		unsigned int MinimumX2(const Eigen::VectorXd &On,
				const Eigen::VectorXd &En,
				Eigen::VectorXd &epsil, double &x2);



		Eigen::VectorXd Jacobian(const Eigen::VectorXd &On,
					 const Eigen::VectorXd &En,
					 const Eigen::VectorXd &epsil);
		Eigen::MatrixXd Hessian(const Eigen::VectorXd &On,
					const Eigen::VectorXd &En,
					const Eigen::VectorXd &epsil);

		void JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil);


		Eigen::MatrixXd Covariance(const Eigen::VectorXd &On,
					   const Eigen::VectorXd &En,
					   const Eigen::VectorXd &epsil);
		Eigen::VectorXd Variance(const Eigen::VectorXd &On,
					 const Eigen::VectorXd &En,
					 const Eigen::VectorXd &epsil);
		
		double X2(const Eigen::VectorXd &On,
			  const Eigen::VectorXd &En,
			  const Eigen::VectorXd &epsil);
		double SysX2(const Eigen::VectorXd &epsil);
		double ObsX2(const Eigen::VectorXd &On,
			     const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);
		Eigen::ArrayXd ObsX2n(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);

		double RawX2(const Eigen::VectorXd &On,
			     const Eigen::VectorXd &En);
		Eigen::ArrayXd RawX2n(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En);

		Eigen::ArrayXd Gamma(const Eigen::VectorXd &En,
				      const Eigen::VectorXd &epsil);
//		Eigen::MatrixXd GammaJac(const Eigen::VectorXd &En,
//				         const Eigen::VectorXd &epsil);
//		Eigen::VectorXd GammaJac(const Eigen::VectorXd &En,
//				         const Eigen::VectorXd &epsil, int j);


		Eigen::ArrayXXd one_F(const Eigen::VectorXd &epsil);
		Eigen::ArrayXd  one_Fk(double err, int k);
		Eigen::ArrayXXd one_Fp(const Eigen::VectorXd &epsil);
		Eigen::ArrayXd  one_Fpk(double err, int k);
		//double F(int k, int n, double eij);
		//double Fp(int k, int n, double eij);
		//double Fp(int k, int n, double dl, double du);


	private:
		//global parameters
		int kVerbosity;
		size_t maxIteration, maxTrials;
		double fitErr;
		bool zeroEpsilons;

		double lm_0, lm_up, lm_down, lm_min;	//control fit parameters

		std::vector<std::shared_ptr<Sample> > _sample;

		int _nBin, _nSys;
		std::map<int, Eigen::ArrayXXd> _sysMatrix;
		Eigen::MatrixXd _corr;
};

#endif
