/*
 * This class defines the matter profile for neutrino propagation
 * it uses Honda prediction for production height
 * and basic trigonometry for neutrino path
 * Earth matter profile is passed via card
 */

#ifndef Atmosphere_H
#define Atmosphere_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

#include <map>
#include <vector>
#include <tuple>

#include <numeric>
#include <algorithm>

#include "tools/CardDealer.h"

#include "physics/Const.h"
#include "physics/Flavours.h"
#include "physics/Oscillator.h"

class Atmosphere
{
	public:
		Atmosphere(CardDealer *cd);

		void LoadProductionHeights();
		double GenerateRandomHeight(Nu::Flavour flv, double cosz, double energy);
		void LoadDensityProfile();
		Oscillator::Profile MatterProfile(double cosz, double atm = 15); // 15 km default

		//std::map<std::string, Eigen::VectorXd>
		//	Oscillate(const std::vector<std::pair<double, double> > &bins, 
		//	const std::map<std::string, std::pair<Nu::Flavour, Nu::Flavour> > &oscf, Oscillator *osc = 0);

	private:
		CardDealer *cd;
		Oscillator::Profile _profile;
		int kVerbosity;

		std::vector<double> _problibs, _energies, _zenithas;
		std::mt19937 gen;

		std::map<int, std::vector<double> > nuE0_table;
		std::map<int, std::vector<double> > nuM0_table;
		std::map<int, std::vector<double> > nuEb_table;
		std::map<int, std::vector<double> > nuMb_table;

};

#endif
