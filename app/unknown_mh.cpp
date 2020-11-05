#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

int main(int argc, char** argv)
{
	if (argc < 4) {
		std::cerr << "Usage: " << argv[0] << " input1 input2 output" << std::endl;
		return 1;
	}

	std::vector<double> x_dCP, tx2_1, cx2_1, tx2_2, cx2_2;
	double dCP, sx2, tx2, cx2;

	std::string line;
	std::ifstream inf(argv[1]);
	while (std::getline(inf, line)) {
		std::stringstream ssl(line);
		ssl >> dCP >> sx2 >> tx2 >> cx2;

		x_dCP.push_back(dCP);
		tx2_1.push_back(tx2);
		cx2_1.push_back(cx2);
	}
	inf.close();

	inf.open(argv[2]);
	while (std::getline(inf, line)) {
		std::stringstream ssl(line);
		ssl >> dCP >> sx2 >> tx2 >> cx2;

		tx2_2.push_back(tx2);
		cx2_2.push_back(cx2);
	}
	inf.close();


	std::ofstream out(argv[3]);
	for (int i = 0; i < x_dCP.size(); ++i)
		out << x_dCP[i] << "\t"
		    << sqrt(std::abs(std::min(tx2_1[i], tx2_2[i])
				   - std::min(cx2_1[i], cx2_2[i]))) << "\n";
	out.close();

	return 0;
}
