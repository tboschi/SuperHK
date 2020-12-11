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
#include <memory>

#include "tools/CardDealer.h"

#include "physics/Const.h"
#include "physics/Flavours.h"
#include "physics/Oscillator.h"

class Atmosphere
{
	public:
		Atmosphere(const std::string &card);
		Atmosphere(const CardDealer &cd);
		Atmosphere(CardDealer *cd);

		void LoadProductionHeights(const CardDealer &cd);
		void LoadDensityProfile(const CardDealer &cd, std::string table_file = "");

		double RandomHeight(Nu::Flavour flv, double cosz, double energy);
		Oscillator::Profile MatterProfile(Nu::Flavour flv, double cosz, double energy);
		Oscillator::Profile MatterProfile(double cosz, double atm = -1);

		//std::map<std::string, Eigen::VectorXd>
		//	Oscillate(const std::vector<std::pair<double, double> > &bins, 
		//	const std::map<std::string, std::pair<Nu::Flavour, Nu::Flavour> > &oscf, Oscillator *osc = 0);

	private:
		Oscillator::Profile _profile;
		int kVerbosity;
		double _atm;

		std::vector<double> _problibs, _energies, _zenithas;
		std::mt19937 gen;

		std::map<int, std::vector<double> > nuE0_table;
		std::map<int, std::vector<double> > nuM0_table;
		std::map<int, std::vector<double> > nuEb_table;
		std::map<int, std::vector<double> > nuMb_table;

};

#endif
