/* Oscillator of the Barger type
 */

#ifndef Oscillator_H
#define Oscillator_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <complex>
#include <utility>

#include "tools/Const.h"
#include "tools/CardDealer.h"
#include "physics/Flavours.h"

#include "Eigen/Dense"

#include "TH1D.h"

class Oscillator
{
	public:
		enum masses
		{
			normal,
			inverted,
			absolute,
		};

		enum pmns
		{
			sin,
			sin2,
			angles,
		};

	public:
		Oscillator(const std::vector<double> &lengths,
			   const std::vector<double> &densities,
			   int numDimension = 3, bool lut = false);
		Oscillator(const std::string &densityFile,
			   int numDimension = 3);
		Oscillator(CardDealer *cd,
			   int numDimension = 3);

		void SetDensity(const std::string &densityFile);

		double Probability(Nu in, Nu out, double energy);
		void Oscillate(Nu in, Nu out, TH1D* h);
		void Reset();
		std::map<double, Eigen::MatrixXd>::iterator FindEnergy(double energy);



		Eigen::MatrixXcd TransitionMatrix(double energy);
		Eigen::MatrixXcd TransitionMatrix(double ff, double l2e);
		void MatterMatrices(Eigen::MatrixXd &dmMatVac,
				    Eigen::MatrixXd &dmMatMat,
				    const double &ff);
		Eigen::VectorXd MatterStates(const double &ff);

		template<masses type>
		void SetMasses(double m1, double m2)
		{
			switch (type)
			{
				case masses::normal:
					SetMasses_NH(m1, m2);
					break;
				case masses::inverted:
					SetMasses_IH(m1, m2);
					break;
				case masses::absolute:
					SetMasses_abs(m1, m2);
					break;
			}
		}
		void SetMasses_NH(double dms21, double dms23);
		void SetMasses_IH(double dms21, double dms23);
		void SetMasses_abs(double ms2, double ms3);

		template<pmns type>
		void SetPMNS(double s12, double s13, double s23, double cp)
		{
			switch (type)
			{
				case pmns::sin:
					SetPMNS_sin(s12, s13, s23, cp);
					break;
				case pmns::sin2:
					SetPMNS_sin2(s12, s13, s23, cp);
					break;
				case pmns::angles:
					SetPMNS_angles(s12, s13, s23, cp);
					break;
			}
		}
		void SetPMNS_sin(double s12, double s13, double s23, double cp);
		void SetPMNS_sin2(double s12, double s13, double s23, double cp);
		void SetPMNS_angles(double t12, double t13, double t23, double cp);

		Eigen::MatrixXcd PMNS();
		Eigen::VectorXd Masses();

	private:
		Eigen::MatrixXcd _pmnsM, pmnsM, transM;
		Eigen::VectorXd dms;
		std::vector<std::pair<double, double> > lens_dens;

		int nDim;	//number of neutrinos
		bool kAntineutrino, kLUT;

		//Fermi constant in SI units times Avogadro's constant (eV² cm³)/(mol GeV)
		double fG = Const::fGF * pow(Const::fhBarC * Const::m/Const::cm, 3) * Const::fNa *
			    Const::GeV/Const::eV;

		std::map<double, Eigen::MatrixXd > mLUT;
		//double fG = 1.52588e-4 / sqrt(8);
};

#endif
