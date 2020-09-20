#include <iostream>

#include "physics/ParameterSpace.h"

int main(int argc, char** argv)
{
	if (argc < 1)
	{
		std::cerr << "Usage: card [parm1, parm2, ...]" << std::endl;
		return 1;
	}

	std::string cardFile = argv[1];
	ParameterSpace *parms = new ParameterSpace(cardFile);

	if (argc < 2)
		std::cout << parms->GetNominalEntry() << std::endl;
	else {
		std::vector<std::string> ps;
		for (int f = 2; f < argc; ++f)
			ps.push_back(std::string(argv[f]));
		std::vector<int> entries = parms->GetScanEntries(ps);
		for (int i : entries)
			std::cout << i << std::endl;
	}

	delete parms;

	return 0;
}
