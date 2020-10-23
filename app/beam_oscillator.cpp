#include <iostream>
#include <fstream>
#include <memory>

#include "tools/CardDealer.h"

#include "physics/Oscillator.h"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Need one parameter: card" << std::endl;
		return 1;
	}

	CardDealer cd(argv[1]);

	std::unique_ptr<Oscillator> osc(new Oscillator(&cd));
	osc->AutoSet(&cd);

	std::string output;
	if (!cd.Get("output", output))
		output = "beam_oscillation.dat";

	if (osc->GetHierarchy() == Oscillator::normal)
		std::cout << "hierarchy is normal\n";
	if (osc->GetHierarchy() == Oscillator::inverted)
		std::cout << "hierarchy is inverted\n";

	std::ofstream out(output.c_str());
	for (double energy = 0; energy <= 10; energy += 0.002) {
		out << energy << "\t";
		out << osc->Probability(Nu::E_, Nu::E_, energy) << "\t";
		out << osc->Probability(Nu::M_, Nu::E_, energy) << "\t";
		out << osc->Probability(Nu::E_, Nu::M_, energy) << "\t";
		out << osc->Probability(Nu::M_, Nu::M_, energy) << "\n";
	}

	out.close();

	return 0;
}
