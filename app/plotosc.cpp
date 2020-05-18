#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

#include "physics/Oscillator.h"
#include "tools/CardDealer.h"
#include "tools/Const.h"

int main(int argc, char** argv)
{
	std::string card = argv[1];

	CardDealer *cd = new CardDealer(card);

	double M12, M23;
	double S12, S13, S23;
	double dCP;

	cd->Get("M12", M12);
	cd->Get("M23", M23);
	cd->Get("S12", S12);
	cd->Get("S13", S13);
	cd->Get("S23", S23);
	cd->Get("dCP", dCP);

	std::vector<double> baseline(1), density(1);

	std::string outFile;
	cd->Get("base", outFile);
	std::ofstream fout(outFile.c_str());

	Oscillator *osc_NH = new Oscillator(baseline, density);

	osc_NH->SetMasses<Oscillator::normal>(M12, M23);
	osc_NH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	double energy = 1.0;
	for (double length = 0.0; length < 5e4; length += 1)
	{
		baseline[0] = length;
		osc_NH->DefineMatter(baseline, density);

		fout << length / energy << "\t"
		     << osc_NH->Probability(Nu::E_, Nu::E_, energy) << "\t"
		     << osc_NH->Probability(Nu::E_, Nu::M_, energy) << "\t"
		     << osc_NH->Probability(Nu::E_, Nu::T_, energy) << "\t"
		     << osc_NH->Probability(Nu::M_, Nu::E_, energy) << "\t"
		     << osc_NH->Probability(Nu::M_, Nu::M_, energy) << "\t"
		     << osc_NH->Probability(Nu::M_, Nu::T_, energy) << std::endl;
	}

	fout.close();

	return 0;
}
