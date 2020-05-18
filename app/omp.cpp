#include <fstream>
#include <iostream>
#include <iomanip>

#include <omp.h>

#include "tools/CardDealer.h"

int main(int argc, char** argv)
{
	std::string cardFile = argv[1];

	CardDealer *cd = new CardDealer(cardFile);

	double fThreads;
	cd->Get("threads", fThreads);
	int nthread = fThreads;
	std::vector<int> vT;

	int size = 1000;

	for (int i = 0; i < size; ++i)
		vT.push_back(i);

#pragma omp parallel for
	for (int i = 0; i < vT.size(); ++i)
	{
		std::cout << "thr " << omp_get_num_threads() << std::endl;
		vT[i] = vT[i] * vT[i];
	}

	std::cout << "vT size " << vT.size() << std::endl;
	for (int i = 0; i < vT.size(); ++i)
		std::cout << i << "\t" << vT.at(i) << std::endl;

	return 0;
}
