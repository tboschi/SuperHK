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

	std::vector<std::vector<double> > allpdf;
	std::vector<std::vector<double> > alldev;
	std::vector<double> sigmas;
	// maybe a loop over scale
	for (double s = -3; s <= 3; s += 0.01) {

		double SKerror = 0.024;
		double scale = 1 + SKerror * s;

		std::vector<double> pdf_scale(nbins, 0.0);
		std::vector<double> derivate(nbins, 0.0);

		// Loop over scaled edges

		for (size_t n = 0; n < nbins - 1; ++n) {

			// unscaled/original bins
			double b0_n = bins[n];
			double b1_n = bins[n + 1];

			//binary search for starting bin
			int m0 = 0;
			int m1 = nbins - 1;
			int mm = (m0 + m1) / 2;
			while (scale * bins[mm] > b0_n
					|| scale * bins[mm+1] < b0_n) {
				if (b0_n <= scale * bins[mm + 1]) {
					m0 = m0;
					m1 = mm;
				}
				else if (b0_n >= scale * bins[mm]) {
					m0 = mm;
					m1 = m1;
				}
				else
					break;
				mm = (m0 + m1) / 2;
			}

			for (int m = mm; m < nbins-1; ++m) {

				// scaled bins
				double b0_m = bins[m];
				double b1_m = bins[m + 1];

				double f = (b1_n - b0_n + scale * (b1_m - b0_m)
					  - std::abs(b0_n - scale * b0_m)
					  - std::abs(b1_n - scale * b1_m) ) / 2.;

				//stop computation because there is no overlap
				if (scale * b0_m > b1_n)
					break;

				double ss = f < 0 ? 0 :
					f > 0 ? 1 :
					0.5;
				//if (f <= 0)
				//{
				//	std::cout << ss << std::endl;
				//	continue;
				//}

				int s0 = b0_n - scale * b0_m < 0 ? -1 : 1;
				int s1 = b1_n - scale * b1_m < 0 ? -1 : 1;
				double fd = (b1_m - b0_m + s0 * b0_m + s1 * b1_m) / 2;

				pdf_scale[n] += pdfs[m] * ss * f / scale
					/ (b1_m - b0_m);

				//derivate[n] += pdfs[m] * SKerror / scale
				//	/ (b1_m - b0_m) * ss * (fd - f / scale);

				derivate[n] += pdfs[m] * pow(SKerror / scale, 2)
					/ (b1_m - b0_m) * 2 * ss * (f / scale - fd);
			}
		}

		if (scale > 1) { //to keep same integral
			double b0_n = bins[nbins - 2];
			double b1_n = bins[nbins - 1];

			// "1-" because this is the bit that is above
			// the maximum (e.g. 30 GeV for 1Rmu)
			// We've already counted the other bit
			double f = 1 - (b1_n - scale * b0_n)
				/ scale / (b1_n - b0_n);
			pdf_scale[nbins - 2] += pdfs[nbins - 2] * f;
		}

		sigmas.push_back(s);
		allpdf.push_back(pdf_scale);
		alldev.push_back(derivate);
	}


	std::ofstream out(argv[2]);

	for (int s = 1; s < sigmas.size() - 1; ++s) {
		out << sigmas[s];
		for (int n = 0; n < nbins; ++n)
			out << "\t" << (allpdf[s+1][n] - 2*allpdf[s][n] + allpdf[s-1][n])
				     / (sigmas[s+1] - sigmas[s]) / (sigmas[s] - sigmas[s-1]);
		out << std::endl;
	}

	out.close();
	out.open(argv[3]);

	for (int s = 0; s < sigmas.size(); ++s) {
		out << sigmas[s];
		for (int n = 0; n < nbins; ++n)
			out << "\t" << alldev[s][n];
		out << std::endl;
	}


	return 0;
}

