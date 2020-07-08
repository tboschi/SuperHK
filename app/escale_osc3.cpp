#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>

#include "TRandom3.h"

int main(int argc, char** argv)
{
	std::ifstream inp(argv[1]);
	std::string line;
	std::vector<double> pdfs, bins;
	while (std::getline(inp, line)) {
		std::stringstream ss(line);
		double bin, cont;
		ss >> bin;
		ss >> cont;
		bins.push_back(bin);
		pdfs.push_back(cont);
	}

	double integral_origin = 0;
	for (double & e : pdfs)
		integral_origin += e;
	std::cout << "original integral " << integral_origin << std::endl;

	int nbins = bins.size();

	// maybe a loop over scale
	double SKerror = 0.024;
	double sigma = std::atof(argv[3]);
	double scale = 1 + SKerror * sigma;
	std::cout << "With sigma " << sigma << " (scale " << scale << ")\n";
	//std::vector<double> bins_scaled = bins;
	//for (double & bin : bins_scaled)
	//	bin *= scale;

	// Calculate the new number of events per bin by seeing
	// how much overlap there is between scaled and unscaled bins.
	// Note: We have underflow and overflow bins here.
	// This is because as we scale the bin edges,
	// events can drop off either end of the scaled binning.
	// These under/overflow bins catch these events,
	// which are then added back on to
	// the first or last bin at the end, respectively.

	std::vector<double> pdf_scale(nbins, 0.0);
	std::vector<double> derivate(nbins, 0.0);
	std::vector<double> hessian(nbins, 0.0);

	// Loop over scaled edges

	for (size_t n = 0; n < nbins - 1; ++n) {

		// unscaled/original bins
		double b0_n = bins[n];
		double b1_n = bins[n + 1];

		//binary search for starting bin
		for (int m = 0; m < nbins-1; ++m) {

			// scaled bins
			double b0_m = bins[m];
			double b1_m = std::min(bins[m + 1], 30 / scale);

			double f = (b1_n - b0_n + scale * (b1_m - b0_m)
					- std::abs(b0_n - scale * b0_m)
					- std::abs(b1_n - scale * b1_m) ) / 2.;

			//stop computation because there is no overlap
			//if (scale * b0_m > b1_n)
			//	break;

			//if (f < 0)
			//	continue;

			double ss = f > 0 ? 1 : 0.5;

			int s0 = b0_n - scale * b0_m < 0 ? -1 : 1;
			int s1 = b1_n - scale * b1_m < 0 ? -1 : 1;
			double fd = (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2;

			f = f > 0 ? f : 0.;
			fd = f > 0 ? fd : 0.;


			pdf_scale[n] += pdfs[m] * ss * f / scale
				/ (b1_m - b0_m);

			derivate[n] += pdfs[m] * SKerror / scale
				/ (b1_m - b0_m) * ss * (fd - f / scale);

			hessian[n] += pdfs[m] * pow(SKerror / scale, 2)
				/ (b1_m - b0_m) * ss * (f / scale - fd);
		}

	}

	//if (scale > 1) { //to keep same integral
	//	double b0_n = bins[nbins - 2];
	//	double b1_n = bins[nbins - 1];
	//	std::cout << "Scaling bin " << b0_n << ":" << b1_n << std::endl;

	//	// "1-" because this is the bit that is above
	//	// the maximum (e.g. 30 GeV for 1Rmu)
	//	// We've already counted the other bit
	//	double f = 1 - (b1_n - scale * b0_n)
	//		/ scale / (b1_n - b0_n);
	//	pdf_scale[nbins - 2] += pdfs[nbins - 2] * f;
	//}


	// Lateral systematics should just migrate events between bins,
	// while keeping the total number of events constant.
	// Check that we haven't changed the total number of events
	double integral_scaled = 0;
	for (double & e : pdf_scale)
		integral_scaled += e;

	//if (std::abs(integral_scaled / integral_origin - 1.) > 1e-6)
	std::cout << std::fixed << std::setprecision(9)
		<< "Normalisation check: "
		<< "original integral = " << integral_origin
		<< ", scaled integral = " << integral_scaled
		<< ", diff = " << integral_origin - integral_scaled
		<< std::endl;

	std::ofstream out(argv[2]);

	for (int n = 0; n < nbins; ++n)
		out << bins[n] << "\t"
			<< pdfs[n] << "\t"
			<< pdf_scale[n] << "\t"
			<< derivate[n] << "\t"
			<< hessian[n] << "\n";

	return 0;
}

