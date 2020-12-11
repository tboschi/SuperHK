#include <iostream>

#include "physics/ParameterSpace.h"

int main(int argc, char** argv)
{
	if (argc < 1)
	{
		std::cerr << "Usage: card [parm1, parm2, ...]" << std::endl;
		return 1;
	}

	CardDealer cd(argv[1]);
	// crate a ParameterSpace object
	// the class flattens a grid of points in a N dimensional space
	// of oscillation parameters, specified in the input card
	std::unique_ptr<ParameterSpace> parms(new ParameterSpace(cd));

	// if only 1 argument is passed, the entry value of the
	// nominal point is shown
	if (argc < 2)
		std::cout << parms->GetNominalEntry() << std::endl;
	// any other argument is a string representing an oscillation parameter
	// the block shows all entries for a scan on the specified parameters
	// when all the others are nominal. Parameters can be
	// M12 M23 S12 S13 S23 CP
	else {
		std::vector<std::string> ps;
		ps.reserve(argc-2);
		for (int f = 2; f < argc; ++f)
			ps.push_back(std::string(argv[f]));
		std::vector<int> entries = parms->GetScanEntries(ps);
		for (int i : entries)
			std::cout << i << std::endl;
	}

	return 0;
}
