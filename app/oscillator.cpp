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

	std::string mh;
	if (!cd->Get("mass_hierarchy", mh))
		mh = "normal";

	Oscillator *osc = new Oscillator(cd);

	if (mh == "normal")
		osc->SetMasses<Oscillator::normal>(M12, M23);
	else if (mh == "inverted")
		osc->SetMasses<Oscillator::inverted>(M12, M23);
	osc->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);

	//osc_NH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);
	//osc_IH->SetPMNS<Oscillator::sin2>(S12, S13, S23, dCP);// * Const::fPi);

	std::string outFile;
	cd->Get("output", outFile);
	std::ofstream fout(outFile.c_str());
	for (double energy = 0.01; energy < 10; energy += 0.01) {
		fout << energy << "\t" << osc->Length()/energy << "\t"
			<< osc->Probability(Nu::M_, Nu::E_, energy) << "\t"
			<< osc->Probability(Nu::M_, Nu::M_, energy) << "\t"
			<< osc->Probability(Nu::M_, Nu::T_, energy) << "\t"
			<< osc->Probability(Nu::Mb, Nu::Eb, energy) << "\t"
			<< osc->Probability(Nu::Mb, Nu::Mb, energy) << "\t"
			<< osc->Probability(Nu::Mb, Nu::Tb, energy) << std::endl;

	}

	fout.close();

	return 0;
}
