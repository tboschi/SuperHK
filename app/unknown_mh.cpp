#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc < 4) {
		std::cerr << "Usage: " << argv[0] << " input1 input2 output" << std::endl;
		return 1;
	}

	std::vector<double> x_dCP, x2_1, x2_2;
	double dCP, sx2;

	std::string line;
	std::ifstream inf(argv[1]);
	while (std::getline(inf, line)) {
		std::stringstream ssl(line);
		ssl >> dCP >> sx2;

		x_dCP.push_back(dCP);
		x2_1.push_back(sx2);
	}
	inf.close();

	inf.open(argv[2]);
	while (std::getline(inf, line)) {
		std::stringstream ssl(line);
		ssl >> dCP >> sx2;

		x2_2.push_back(sx2);
	}
	inf.close();


	std::ofstream out(argv[3]);
	for (int i = 0; i < x_dCP.size(); ++i)
		out << x_dCP[i] << "\t"
		    << std::min(x2_1[i], x2_2[i]) << "\n";
	out.close();

	return 0;
}
