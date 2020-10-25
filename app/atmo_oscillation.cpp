#include <iostream>
#include <fstream>

//#include "physics/Oscillator.h"
#include <iostream>
#include <fstream>
#include <memory>

#include "tools/CardDealer.h"

#include "physics/Atmosphere.h"
#include "physics/Oscillator.h"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Need one parameter: card" << std::endl;
		return 1;
	}

	CardDealer cd(argv[1]);

	std::unique_ptr<Atmosphere> atm(new Atmosphere(&cd));
	std::unique_ptr<Oscillator> osc(new Oscillator(&cd));
	osc->AutoSet(&cd);

	std::string output;
	if (!cd.Get("output", output))
		output = "atmospheric_oscillation.dat";
	std::ofstream out(output.c_str());
	for (double cosz = -1; cosz <= 1; cosz += 0.002) {
		Oscillator::Profile pp = atm->MatterProfile(cosz);
		osc->SetMatterProfile(pp);
		for (double lener = -1; lener <= log10(300); lener += 0.002) {
			double energy = pow(10, lener);
			out << cosz << "\t" << energy << "\t";
			out << osc->Probability(Nu::E_, Nu::E_, energy) << "\t";
			out << osc->Probability(Nu::M_, Nu::E_, energy) << "\t";
			out << osc->Probability(Nu::E_, Nu::M_, energy) << "\t";
			out << osc->Probability(Nu::M_, Nu::M_, energy) << "\n";
		}
	}

	out.close();

	return 0;
}
