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
#include <tuple>
#include <algorithm>
#include <numeric>

#include "tools/CardDealer.h"
#include "physics/Const.h"
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

		typedef std::vector<std::tuple<double, double, double> > Profile;
		// anyone can access this without object
		static Oscillator::Profile GetMatterProfile(const std::string &densityFile);

		static double Length(const Oscillator::Profile &ld);	// return total baseline
		static double Density(const Oscillator::Profile &ld);	// return average density
		static double ElectronDensity(const Oscillator::Profile &ld); // return average electron density


		double Length();	// return total baseline
		double Density();	// return average density
		double ElectronDensity(); // return average electron density

		void SetMatterProfile(const Oscillator::Profile &l_d);


	public:
		Oscillator(const std::vector<double> &lengths,
			   const std::vector<double> &densities,
			   bool lut = false, double threshold = 1e-9);
		Oscillator(const std::vector<double> &lengths,
			   const std::vector<double> &densities,
			   const std::vector<double> &electrons,
			   bool lut = false, double threshold = 1e-9);

		Oscillator(const std::string &densityFile,
			   bool lut, double threshold);
		Oscillator(const std::string &cd);
		Oscillator(CardDealer *cd);

		void FromCard(CardDealer *cd);

		double Probability(Nu::Flavour in, Nu::Flavour out,
				double energy, bool force = false);
		//void Oscillate(Nu in, Nu out, TH1D* h);
		Eigen::VectorXd Oscillate(Nu::Flavour in, Nu::Flavour out,
				const std::vector<double> &bins);
		void Reset();
		std::map<double, Eigen::MatrixXd>::iterator FindEnergy(double energy);

		Eigen::MatrixXcd TransitionMatrix(double energy);
		Eigen::MatrixXcd TransitionMatrix(double ff, double l2e);
		void MatterMatrices(Eigen::MatrixXd &dmMatVac,
				    Eigen::MatrixXd &dmMatMat,
				    double ff);
		Eigen::VectorXd MatterStates(double ff, int off = 0);

		void AutoSet(CardDealer *cd);

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

			//reset lookup table
			Reset();
		}
		void SetMasses_NH(double dms21, double dms23);
		void SetMasses_IH(double dms21, double dms23);
		void SetMasses_abs(double ms2, double ms3);
		masses GetHierarchy();

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

			//reset lookup table
			Reset();
		}
		void SetPMNS_sin(double s12, double s13, double s23, double cp);
		void SetPMNS_sin2(double s12, double s13, double s23, double cp);
		void SetPMNS_angles(double t12, double t13, double t23, double cp);

		Eigen::MatrixXcd PMNS();
		Eigen::VectorXd Masses();

	private:
		Eigen::MatrixXcd _pmns; //, _pmnsM, pmnsM, trans;
		Eigen::VectorXd dms;
		Profile _lens_dens;

		int _dim;	//number of neutrinos
		double _thr;
		bool kLUT;

		//Fermi constant in SI units times Avogadro's constant (eV² cm³)/(mol GeV)
		//double fG2 = Const::GF * pow(Const::hBarC * Const::m/Const::cm, 3) * Const::Na * Const::GeV/Const::eV;
		//double fG = Const::GF * Const::Na * pow(Const::hBarC, 3) * 1e13;
		double fG, fG2;

		std::map<double, Eigen::MatrixXd > mLUT;
		//double fG = 1.52588e-4 / sqrt(8);
};

#endif
