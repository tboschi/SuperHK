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
#include <memory>

#include "tools/CardDealer.h"

class ParameterSpace
{
	public:

		using Binning = std::map<std::string, std::vector<double> >;

		ParameterSpace(const std::string &card);
		ParameterSpace(const CardDealer &cd);
		ParameterSpace(CardDealer *cd);

		void Init(const Binning &parms);

		Binning GetHistogramBinning();
		Binning GetBinning();

		int GetEntries();
		double GetPenalty(int n);
		void GetEntry(int n, double &M12, double &M23,
			      double &S12, double &S13, double &S23, double &dCP);
		std::map<std::string, double> GetEntry(int n);

		void GetNominal(double &M12, double &M23,
			      double &S12, double &S13, double &S23, double &dCP);
		std::map<std::string, double> GetNominal();
		int GetNominalEntry();
		std::vector<int> GetScanEntries(const std::vector<std::string> &p);

	private:
		//std::map<std::string, double*> varmap;	//map to address of variables
		//std::map<std::string, double*>::iterator iv;
		Binning _binning, _penals;	//map of bins
		std::map<std::string, int> _nominal;
};

#endif

