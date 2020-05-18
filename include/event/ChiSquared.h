/* ChiSquared
 * compute the likelihood ratio chi-squared using the pull approach to include systematics
 */

#ifndef ChiSquared_H
#define ChiSquared_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixT.h"
#include "TKey.h"

#include "Eigen/Dense"
#include "Eigen/SVD"

#include "tools/CardDealer.h"
#include "event/Reco.h"
#include "physics/Oscillator.h"

class ChiSquared
{
	public:
		ChiSquared(CardDealer *card, int nbins = -1);

		void Init();
		void LoadCorrelation();
		void LoadSystematics();

		Eigen::VectorXd ConstructSpectrum(Oscillator *osc);
		Eigen::VectorXd LoadSpectrum(int pt, std::string from = "");

		bool PointInFile(TFile *f, int pt);

		int NumSys();
		int NumBin();
		int DOF();

		Eigen::VectorXd FitX2(const Eigen::VectorXd &On,
				      const Eigen::VectorXd &En);
		bool FindMinimum(const Eigen::VectorXd &On,
				 const Eigen::VectorXd &En,
				 Eigen::VectorXd &epsil);

		Eigen::VectorXd Jacobian(const Eigen::VectorXd &On,
					 const Eigen::VectorXd &En,
					 const Eigen::VectorXd &epsil);
		Eigen::MatrixXd Hessian(const Eigen::VectorXd &On,
					const Eigen::VectorXd &En,
					const Eigen::VectorXd &epsil);
		void JacobianHessian(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil);

		void JacobianHessian2(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
				     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
				     const Eigen::VectorXd &epsil);

		Eigen::MatrixXd Covariance(const Eigen::VectorXd &On,
					   const Eigen::VectorXd &En,
					   const Eigen::VectorXd &epsil);
		
		double X2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			  const Eigen::VectorXd &epsil);
		double ObsX2(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);
		std::vector<double> ObsX2n(const Eigen::VectorXd &On, const Eigen::VectorXd &En,
			     const Eigen::VectorXd &epsil);
		double SysX2(const Eigen::VectorXd &epsil);

		Eigen::VectorXd Gamma(const Eigen::VectorXd &En,
				      const Eigen::VectorXd &epsil);
		Eigen::MatrixXd GammaJac(const Eigen::VectorXd &En,
				         const Eigen::VectorXd &epsil);
		Eigen::VectorXd GammaJac(const Eigen::VectorXd &En,
				         const Eigen::VectorXd &epsil, int j);
		Eigen::VectorXd GammaHes(const Eigen::VectorXd &En,
				         const Eigen::VectorXd &epsil, int j, int i);


		double F(int k, int n, double eij);
		double Fp(int k, int n, double eij);
		double Fp(int k, int n, double dl, double du);

		int StartingBin(int n, double scale);

		//function to compute the factor
		typedef double (ChiSquared::*FactorFn)(const std::vector<double> &);

		double Scale(FactorFn factor,
			     const Eigen::VectorXd &En, double sigma, int t);
		double Scale(const Eigen::VectorXd &En, double sigma, int t);
		double ScaleJac(const Eigen::VectorXd &En, double sigma, int t);
		double ScaleHes(const Eigen::VectorXd &En, double sigma, int t);

		double ScaleNor(const std::vector<double> &term);
		double ScaleJac(const std::vector<double> &term);
		double ScaleHes(const std::vector<double> &term);




	private:
		CardDealer *cd;

		//global parameters
		int kVerbosity;
		int badFitThreshold;
		int maxIteration;
		double err;

		// for loops
		unsigned int kMode;
		unsigned int kChan;
		unsigned int kType;
		unsigned int kHorn;
		unsigned int kOscf;
		unsigned int kFIn ;
		unsigned int kFOut;

		std::vector<std::string> mode;
		std::vector<std::string> chan;
		std::vector<std::string> type;
		std::vector<std::string> horn;
		std::vector<std::string> oscf;
		std::vector<Nu> fin;
		std::vector<Nu> fout;

		std::map<std::string, double> scale;
		std::map<std::string, double>::iterator is;
		unsigned int kScale;

		///pulls
		//Eigen::VectorXd epsil;
		//systematics, corrlation matrix and fij
		Eigen::MatrixXd corr;

		//promote all systematics to be spline
		std::map<int, Eigen::MatrixXd> sysMatrix;

		//number of bins
		int _nBin;
		std::vector<double> _bins;

		TFile *spectrumFile;
		std::string input;
};

#endif
