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

	std::string densityFile;
	cd->Get("density_profile", densityFile);

	std::string outFile;
	cd->Get("base", outFile);
	std::ofstream fout(outFile.c_str());

	Oscillator *osc_NH = new Oscillator(densityFile);
	Oscillator *osc_IH = new Oscillator(densityFile);

	osc_NH->SetMasses<Oscillator::normal>(M12, M23);
	osc_IH->SetMasses<Oscillator::inverted>(M12, M23);

	//osc_NH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);
	//osc_IH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);

	double energy = 0.6;

	for (double S23 = 0.4; S23 < 0.65; S23 += 0.1)
	{
		//for (double dCP = -Const::fPi; dCP < Const::fPi+1e-9; dCP += 2*Const::fPi / 500)
		for (double dCP = -Const::fPi; dCP < Const::fPi+1e-9; dCP += Const::fPi / 2. )
		{
			osc_NH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);
			osc_IH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

			fout << dCP << "\t" << S23 << "\t"
				<< osc_NH->Probability(Nu::M_, Nu::E_, energy) << "\t"
				<< osc_NH->Probability(Nu::Mb, Nu::Eb, energy) << "\t"
				<< osc_IH->Probability(Nu::M_, Nu::E_, energy) << "\t"
				<< osc_IH->Probability(Nu::Mb, Nu::Eb, energy) << std::endl;

		}

		fout << std::endl << std::endl;
	}

	fout.close();

	return 0;
}
