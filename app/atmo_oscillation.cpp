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

	// create a CardDealer object
	// this will load the parameters
	CardDealer cd(argv[1]);

	// an Atmosphere object can be quite heavy
	// it is best to create it dynamically on the head
	std::unique_ptr<Atmosphere> atm(new Atmosphere(cd));

	// see comments on app/beam_oscillaiton.cpp
	std::unique_ptr<Oscillator> osc(new Oscillator(cd));
	osc->AutoSet(cd);

	std::string output;
	if (!cd.Get("output", output))
		output = "atmospheric_oscillation.dat";
	std::ofstream out(output.c_str());

	// nested loop over cosine of Zenith angle and
	// neutrino energy. from different angles the
	// neutrino sees different Earth density profiles
	for (double cosz = -1; cosz <= 1; cosz += 0.002) {
		// create a profile at this angle and
		// at height specified in card
		Oscillator::Profile pp = atm->MatterProfile(cosz);
		osc->SetMatterProfile(pp);
		for (double lener = -1; lener <= log10(300); lener += 0.002) {
			double energy = pow(10, lener);
			// production height can also be generated from
			// the Honda model, but it depends on the neutrino
			// flavour and its energy
			// Oscillator::Profile pp = atm->MatterProfile(Nu::E_, cosz, energy);
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
