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
	TRandom3 *rnd = new TRandom3(0);
	std::ifstream inp(argv[1]);
	std::string line;
	std::vector<double> pdf_origin, bins_origin;
	while (std::getline(inp, line)) {
		std::stringstream ss(line);
		double bin, cont;
		ss >> bin;
		ss >> cont;
		bins_origin.push_back(bin);
		//pdf_origin.push_back(cont*(0.75 + 0.5*rnd->Rndm()));
		pdf_origin.push_back(cont);
	}

	pdf_origin.pop_back();
	double integral_origin = 0;
	for (double & pdf : pdf_origin)
		integral_origin += pdf;

	int nbins = bins_origin.size();

	// maybe a loop over scale
	
	std::vector<std::vector<double> > allpdf;
	std::vector<double> sigmas;

	for (double sigma = -3; sigma <= 3; sigma += 1) {
		double scale = 1 + 0.024 * sigma;
		std::cout << "With sigma " << sigma
			  << " (scale " << scale << ")\n";
		std::vector<double> bins_scaled = bins_origin;
		for (double & bin : bins_scaled)
			bin *= scale;

		// Calculate the new number of events per bin by seeing
		// how much overlap there is between scaled and unscaled bins.
		// Note: We have underflow and overflow bins here.
		// This is because as we scale the bin edges,
		// events can drop off either end of the scaled binning.
		// These under/overflow bins catch these events,
		// which are then added back on to
		// the first or last bin at the end, respectively.

		std::vector<double> pdf_scaled(nbins, 0.0);
		// Loop over scaled edges
		for (size_t s = 0; s < nbins - 1; ++s) {

			// scaled bins
			double bin_lo_s = bins_scaled[s];
			double bin_up_s = bins_scaled[s + 1];

			// Loop over unscaled/original edges
			for (size_t o = 0; o < nbins - 1; ++o) {

				// unscaled/original bins
				double bin_lo_o = bins_origin[o];
				double bin_up_o = bins_origin[o + 1];

				double factor = 0.0;

				// Unscaled bin entirely larger than scaled bin,
				// stop checking higher unscaled
				if (bin_lo_o > bin_up_s)
					break;
				// Scaled bin entirely larger than unscaled bin,
				// skip to next unscaled bin
				else if (bin_up_o < bin_lo_s)
					continue;
				// Unscaled bin subset of scaled bin
				else if (bin_lo_o >= bin_lo_s && bin_up_o <= bin_up_s)
					factor = (bin_up_o - bin_lo_o)
					       / (bin_up_s - bin_lo_s);
				// Scaled bin subset of unscaled bin
				// Take the contents of this bin unchanged,
				// there could be additional scaled bins contributing
				else if (bin_lo_o <= bin_lo_s && bin_up_o >= bin_up_s)
					factor = 1.0;
				// Scaled bin overlaps low part of unscaled bin
				else if (bin_up_o > bin_up_s && bin_lo_o >= bin_lo_s)
					factor = (bin_up_s - bin_lo_o)
					       / (bin_up_s - bin_lo_s);
				// Scaled bin overlaps high part of unscaled bin
				else if (bin_lo_o < bin_lo_s && bin_up_o <= bin_up_s)
					factor = (bin_up_o - bin_lo_s)
					       / (bin_up_s - bin_lo_s);
				// Shouldn't get here
				else {
					std::cout << "here\n";
					factor = 0.0;
				}

				if (factor <= 0.0) {
					std::cerr << "Got a factor <= 0 for " << s
						  << ", " << o
						  << ".This shouldn't happen!"
						  << std::endl;
					continue;
				}

				//inverted bins?
				pdf_scaled[o] += pdf_origin[s] * factor;
			}
		}

		// add the overflow
		// don't add an underflow,
		// because the first couple of bins are actually empty
		if (scale > 1) {
			double bin_lo_s = bins_scaled[nbins - 2];
			double bin_up_s = bins_scaled[nbins - 1];

			double bin_lo_o = bins_origin[nbins - 2];
			double bin_up_o = bins_origin[nbins - 1];

			// "1-" because this is the bit that is above
			// the maximum (e.g. 30 GeV for 1Rmu)
			// We've already counted the other bit
			double factor = 1 - (bin_up_o - bin_lo_s)
					  / (bin_up_s - bin_lo_s);
			pdf_scaled[nbins - 2] += pdf_origin[nbins - 2] * factor;
		}

		// Lateral systematics should just migrate events between bins,
		// while keeping the total number of events constant.
		// Check that we haven't changed the total number of events
		double integral_scaled = 0;
		for (double & pdf : pdf_scaled)
			integral_scaled += pdf;

		//if (std::abs(integral_scaled / integral_origin - 1.) > 1e-6)
			std::cout << std::fixed << std::setprecision(9)
				  << "Normalisation check: "
				  << "original integral = " << integral_origin
				  << ", scaled integral = " << integral_scaled
				  << ", diff = " << integral_origin - integral_scaled
				  << std::endl;

		sigmas.push_back(sigma);
		allpdf.push_back(pdf_scaled);
	}

	std::ofstream out(argv[2]);

	for (int n = 0; n < nbins; ++n) {
		out << bins_origin[n] << "\t" << pdf_origin[n];
		for (int s = 0; s < allpdf.size(); ++s)
			out << "\t" << allpdf[s][n];
		out << "\n";
	}
	return 0;
}
