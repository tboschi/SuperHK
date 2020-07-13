/* ParameterSpace
 * builds the parameter multidimensional space
 */

#ifndef ParameterSpace_H
#define ParameterSpace_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iterator>
#include <cmath>

#include "tools/CardDealer.h"

class ParameterSpace
{
	public:

		typedef std::map<std::string, std::vector<double> > Binning;

		ParameterSpace(CardDealer *card);

		void Init();

		Binning GetHistogramBinning();
		Binning GetBinning();

		int GetEntries();
		double GetPenalty(int n);
		void GetEntry(int n, double &M12, double &M23,
			      double &S12, double &S13, double &S23, double &dCP);
		std::map<std::string, double> GetEntry(int n);

	private:
		CardDealer *cd;

		//std::map<std::string, double*> varmap;	//map to address of variables
		//std::map<std::string, double*>::iterator iv;
		Binning binning, penals;	//map of bins
};

#endif

