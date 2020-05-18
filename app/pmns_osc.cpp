#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	//set up oscillation
	std::string densityFile;
	cd->Get("density_profile", densityFile);

	double M12, M23;
	double S12, S13, S23;
	double dCP;

	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	double baseline;	//km
	double density;		//g/cm³
	cd->Get("baseline", baseline);
	cd->Get("density",  density);

	double kanti;
	cd->Get("anti", kanti);
	bool anti = static_cast<bool>(kanti);

	std::string oscillator, prob3, matrix;
	cd->Get("output1", oscillator);
	cd->Get("output2", prob3);
	cd->Get("matrix", matrix);
	std::ofstream outOsc(oscillator.c_str());
	std::ofstream mat(matrix.c_str());

	std::vector<double> vL(1, baseline), vD(1, density);
	Oscillator *osc = new Oscillator(vL, vD);
	osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);


	Nu fin_[3]  = {Nu::E_, Nu::M_, Nu::T_};
	Nu fout_[3] = {Nu::E_, Nu::M_, Nu::T_};
	Nu finb[3]  = {Nu::Eb, Nu::Mb, Nu::Tb};
	Nu foutb[3] = {Nu::Eb, Nu::Mb, Nu::Tb};

	for (double en = 0.00; en < 5.0; en += 0.01)
	{

		outOsc << en << "\t";
		mat << std::setfill(' ') << std::setprecision(5) << std::fixed;
		mat << std::setw(8) << en << "\n";

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (anti)
				{
					mat << osc->Probability(finb[i], foutb[j], en) << "\t";
					outOsc << "\t" << osc->Probability(finb[i], foutb[j], en);
				}
				else
				{
					mat << osc->Probability(fin_[i], fout_[j], en) << "\t";
					outOsc << "\t" << osc->Probability(fin_[i], fout_[j], en);
				}
			}
		}
		mat << "\n\n";
		outOsc << "\n";
	}
	//init_mixing_matrix(m21, lm32, sqrt(S12), sqrt(S23), sqrt(S13), dCP);
	/*
	std::cout << pow(mix[0][0][0], 2) + pow(mix[0][0][1], 2) << "\t"
		  << pow(mix[0][1][0], 2) + pow(mix[0][1][1], 2) << "\t"
		  << pow(mix[0][2][0], 2) + pow(mix[0][2][1], 2) << std::endl;
	std::cout << pow(mix[1][0][0], 2) + pow(mix[1][0][1], 2) << "\t"
		  << pow(mix[1][1][0], 2) + pow(mix[1][1][1], 2) << "\t"
		  << pow(mix[1][2][0], 2) + pow(mix[1][2][1], 2) << std::endl;
	std::cout << pow(mix[2][0][0], 2) + pow(mix[2][0][1], 2) << "\t"
		  << pow(mix[2][1][0], 2) + pow(mix[2][1][1], 2) << "\t"
		  << pow(mix[2][2][0], 2) + pow(mix[2][2][1], 2) << std::endl;
		  */

	return 0;
}
