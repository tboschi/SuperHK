#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "tools/CardDealer.h"

int main(int argc, char** argv)
{
	CardDealer cd(argv[1]);

	std::vector<std::string> allkeys = cd.ListKeys();

	bool collectX = true;
	std::vector<double> xaxis;
	std::vector<std::vector<double> > ylines;
	std::vector<std::string> names;
	for (const std::string &k : allkeys) {
		std::string file;
		if (!cd.Get(k, file))
			continue;
	
		std::vector<double> yaxis;

		std::ifstream in(file.c_str());
		std::string line;
		while (std::getline(in, line))
		{
			if (line.find_first_of('#') != std::string::npos)
				line.erase(line.find_first_of('#'));

			if (line.empty())
				continue;

			double word;
			std::vector<double> words;
			std::stringstream ssl(line);
			while (ssl >> word)
				words.push_back(word);

			if (collectX && words.size())	
				xaxis.push_back(words.front());

			if (words.size() > 1)
				yaxis.push_back(words.at(1));
		}

		collectX = !xaxis.size();	//if xaxis is full, stop collecting x
		ylines.push_back(yaxis);
		names.push_back(k);
	}

	std::ofstream aout("exclusion_all.dat");
	std::ofstream dout("exclusion_diff.dat");
	aout << "#";
	dout << "#";
	bool add = false;
	for (const std::string &k : names) {
		aout << "\t" << k;
		if (add)
			dout << "\t" << k;
		add = true;
	}
	aout << std::endl;
	dout << std::endl;

	for (size_t i = 0; i < xaxis.size(); ++i)
	{
		aout << xaxis[i];
		dout << xaxis[i];

		bool add = false;
		for (const auto &yl : ylines) {
			aout << "\t" << yl[i];
			if (add)
				dout << "\t" << ylines.front()[i] - yl[i];
			add = true;
		}

		aout << std::endl;
		dout << std::endl;
	}

	aout.close();
	dout.close();

	return 0;
}
