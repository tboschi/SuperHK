#include <iostream>
#include <getopt.h>

#include "tools/CardDealer.h"

void Usage(char* name);

int main(int argc, char** argv)
{
	CardDealer *cd = new CardDealer(argv[1], true);
	std::vector<std::string> string_keys = cd->ListStringKeys();
	std::vector<std::string> double_keys = cd->ListDoubleKeys();

	for (int i = 0; i < string_keys.size(); ++i)
	{
		std::vector<std::string> vs;
		std::cout << string_keys.at(i) << "\n";
		if (cd->Get(string_keys.at(i), vs))
			for (int j = 0; j < vs.size(); ++j)
				std::cout << "\t" << vs[j] << std::endl;
		std::cout << std::endl;
	}

	for (int i = 0; i < double_keys.size(); ++i)
	{
		std::vector<double> vd;
		std::cout << double_keys.at(i) << "\n";
		if (cd->Get(double_keys.at(i), vd))
			for (int j = 0; j < vd.size(); ++j)
				std::cout << "\t" << vd[j] << std::endl;
		std::cout << std::endl;
	}

	return 0;
}
