#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

int main(int argc, char** argv)
{
	std::vector<double> axis, sens1, sens2;

	std::ifstream in(argv[1]);
	std::string line;
	while (std::getline(in, line))
	{
		std::stringstream ssl(line);
		std::vector<double> words;
		double word;

		while (ssl >> word)
			words.push_back(word);

		axis.push_back(words.at(0));
		sens1.push_back(words.at(1));
	}

	in.close();

	in.open(argv[2]);
	while (std::getline(in, line))
	{
		std::stringstream ssl(line);
		std::vector<double> words;
		double word;

		while (ssl >> word)
			words.push_back(word);

		sens2.push_back(std::sqrt(words.at(1)));
	}

	std::ofstream out(argv[3]);
	for (int i = 0; i < axis.size(); ++i)
		out << axis[i] << "\t" << std::min(sens1[i], sens2[i]) << std::endl;

	out.close();

	return 0;
}
