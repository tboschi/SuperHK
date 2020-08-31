#include <fstream>
#include <iostream>

#include "physics/ParameterSpace.h"

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cerr << "Need two parameters: card output" << std::endl;
		return 1;
	}

	std::string cardFile = argv[1], output = argv[2];
	ParameterSpace *parms = new ParameterSpace(cardFile);

	std::ofstream out(output);
	out << parms->GetNominalEntry() << std::endl;
	out.close();

	delete parms;

	return 0;
}
