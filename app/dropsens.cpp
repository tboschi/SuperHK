#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

int main(int argc, char** argv)
{
	std::cout << "Reading " << argc-1 << " files" << std::endl;

	bool collectX = true;
	std::vector<double> xaxis;
	std::map<int, std::vector<double> > ylines;
	for (int f = 1; f < argc; ++f)
	{
		std::vector<double> yaxis;

		std::ifstream in(argv[f]);
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
		ylines[f] = yaxis;
	}

	std::ofstream aout("exclusion_all.dat");
	std::ofstream dout("exclusion_diff.dat");
	aout << "#";
	dout << "#";
	for (int f = 1; f < argc; ++f)
	{
		std::string name(argv[f]);
		name.erase(0, name.find("errorstudy/")+11);
		name.erase(name.find_first_of('/'));
		aout << "\t" << name;
		if (f > 1)
			dout << "\t" << name;
	}
	aout << std::endl;
	dout << std::endl;

	for (int i = 0; i < xaxis.size(); ++i)
	{
		aout << xaxis[i];
		dout << xaxis[i];

		for (int f = 1; f < argc; ++f)
		{
			aout << "\t" << ylines[f][i];
			if (f > 1)
				dout << "\t" << ylines[1][i] - ylines[f][i];
		}

		aout << std::endl;
		dout << std::endl;
	}

	aout.close();
	dout.close();

	return 0;
}
