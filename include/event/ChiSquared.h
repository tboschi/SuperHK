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
#include <utility>
#include <numeric>

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
		void DefineBinning();
		void LoadCorrelation();
		void LoadSystematics();

		Eigen::VectorXd ConstructSpectrum(Oscillator *osc = 0);
		Eigen::VectorXd LoadSpectrum(int pt, std::string from = "");
		std::map<std::string, Eigen::MatrixXd> BuildSpectrum(Oscillator *osc = 0);

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

		//void JacobianHessian2(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
		//		     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		//		     const Eigen::VectorXd &epsil);
		//void JacobianHessian3(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
		//		     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		//		     const Eigen::VectorXd &epsil);
		//void JacobianHessian4(Eigen::VectorXd &jac, Eigen::MatrixXd &hes,
		//		     const Eigen::VectorXd &On, const Eigen::VectorXd &En,
		//		     const Eigen::VectorXd &epsil);

		Eigen::MatrixXd Covariance(const Eigen::VectorXd &On,
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

		Eigen::ArrayXd Gamma(const Eigen::VectorXd &En,
				      const Eigen::VectorXd &epsil);
//		Eigen::MatrixXd GammaJac(const Eigen::VectorXd &En,
//				         const Eigen::VectorXd &epsil);
//		Eigen::VectorXd GammaJac(const Eigen::VectorXd &En,
//				         const Eigen::VectorXd &epsil, int j);


		Eigen::ArrayXXd one_F(const Eigen::VectorXd &epsil);
		Eigen::ArrayXd one_Fk(double err, int k);
		Eigen::ArrayXXd one_Fp(const Eigen::VectorXd &epsil);
		Eigen::ArrayXd one_Fpk(double err, int k);
		//double F(int k, int n, double eij);
		//double Fp(int k, int n, double eij);
		//double Fp(int k, int n, double dl, double du);

		int StartingBin(std::string it, double shift, int n);
		int EndingBin(std::string it, double shift, int n);
		std::string TypeFromBin(int n);

		//function to compute the factor
		typedef double (ChiSquared::*FactorFn)(const std::vector<double> &);
		typedef double (ChiSquared::*LazyEvent)(int n);

		//double GammaHes(int n, int k = -1, double ek = 0, int j = -1, double ej = 0);

		std::vector<std::pair<int, int> >
			AllSlices(std::string it, double skerr);
		std::vector<Eigen::ArrayXd> AllScale(FactorFn factor,
						 std::string it, double skerr);

		//double QuickScale(FactorFn factor, const Eigen::VectorXd &En,
		//		//const Eigen::VectorXd &gh,
		//		const std::vector<double> &global,
		//		double b0_n, double b1_n,
		//		double scale_err, double shift,
		//		int m0, int m1, int off);
		//		  //int k = -1, double ek = 0, int j = -1, double ej = 0);

		double Scale(FactorFn factor,
			     const Eigen::ArrayXd &En,
			     double skerr, int n, std::string it = "", int m0 = -1);

		double Scale(const Eigen::ArrayXd &En, 
			     double skerr, int n, std::string it = "", int m0 = -1);
		double Scale(const Eigen::VectorXd &En, 
			     double skerr, int n, std::string it = "", int m0 = -1);

			     //int k = -1, double ek = 0, int j = -1, double = 0);
		//double ScaleJac(const Eigen::VectorXd &En, 
		//		double skerr, int n, std::string it = "", int m0 = -1);
		//		//int k = -1, double ek = 0, int j = -1, double = 0);
		//double ScaleHes(const Eigen::VectorXd &En, 
		//		const Eigen::VectorXd &sk,
		//		int n, std::string it = "", int m0 = -1);
		//		//int k = -1, double ek = 0, int j = -1, double = 0);

		double Nor(const std::vector<double> &term);
		double Jac(const std::vector<double> &term);
		double Hes(const std::vector<double> &term);




	private:
		CardDealer *cd;

		//global parameters
		int kVerbosity;
		int badFitThreshold;
		int maxIteration;
		double fitErr;
		bool zeroEpsilons;

		double lm_0, lm_up, lm_down;	//control fit parameters

		std::vector<std::string> _mode;
		std::vector<std::string> _chan;
		std::vector<std::string> _type;
		std::vector<std::string> _horn;
		std::vector<std::string> _oscf;
		std::vector<Nu> _fin;
		std::vector<Nu> _fout;

		std::map<std::string, double> _scale;
		std::map<std::string, int> _type_scale;
		unsigned int _nScale;

		///pulls
		//Eigen::VectorXd epsil;
		//systematics, corrlation matrix and fij
		Eigen::MatrixXd corr;

		//promote all systematics to be spline
		std::map<int, Eigen::ArrayXXd> sysMatrix;

		//number of bins
		int _nBin, _allBin, _nSys;
		std::map<std::string, std::pair<int, int> > _binpos;
		std::map<std::string, std::pair<int, int> > _limits;
		std::map<std::string, std::vector<double> > _global;
		// global binning from Eigen binning
		std::string _tlim;

		TFile *spectrumFile;
		std::string input;
};

#endif
