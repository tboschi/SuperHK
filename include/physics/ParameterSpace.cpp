#include "ParameterSpace.h"

ParameterSpace::ParameterSpace(const std::string &card)
{
	CardDealer cd(card);

	Binning parms;
	if (!cd.Get("parm_", parms))
		throw std::invalid_argument("ParameterSpace: no parameters in oscillation card\n");
	cd.Get("penalty_", _penals);

	Init(parms);
}

ParameterSpace::ParameterSpace(const CardDealer &cd)
{
	Binning parms;
	if (!cd.Get("parm_", parms))
		throw std::invalid_argument("ParameterSpace: no parameters in oscillation card\n");
	cd.Get("penalty_", _penals);

	Init(parms);
}

ParameterSpace::ParameterSpace(CardDealer *cd)
{
	Binning parms;
	cd->Get("parm_", parms);
	cd->Get("penalty_", _penals);

	Init(parms);
}

void ParameterSpace::Init(const Binning &parms)
{
	//each parameters has start, end and points saved in a vector with three elems
	//bin width = bw = (e - s) / (p - 1)
	//Binning::iterator ip;	//loop on parms
	//for (ip = parms.begin(); ip != parms.end(); ++ip)
	for (const auto &ip : parms)
	{
		//if less than three elements skip 
		//if (ip->second.size() < 3)
		if (ip.second.size() < 3)
		{
			std::vector<double> bins(1);
			bins[0] = ip.second[0];

			//varmap[ip->first] = new double;
			_binning[ip.first] = bins;
			_nominal[ip.first] = 0;
		}
		else
		{
			if (ip.second.size() > 3)	//nominal point
				_nominal[ip.first] = int(ip.second[3]);
			else				// take middle point
				_nominal[ip.first] = int(ip.second[2]-1) / 2;

			//create array with bins
			std::vector<double> bins(ip.second[2]);

			//bin width
			double bw = std::abs(ip.second[1] - ip.second[0]) / 
				    (ip.second[2] - 1);

			//fill first bin
			bins[0] = std::min(ip.second[0], ip.second[1]);

			//and add step until end
			for (size_t i = 1; i < ip.second[2]; ++i)
				bins[i] = bins[i-1] + bw;

			//varmap[ip->first] = new double;		//null pointer
			_binning[ip.first] = bins;
		}
	}
}

ParameterSpace::Binning ParameterSpace::GetBinning()
{
	return _binning;
}

//the binning should be centered, so given s(tart), e(nd), and p(oints)
//we want s - bw/2, s + bw/2, s + 3 bw/2, ... e - bw /2, e + bw/2
//bin width = bw = (e - s) / (p - 1) and number of bins is p
ParameterSpace::Binning ParameterSpace::GetHistogramBinning()
{
	Binning histbin;
	//for (ib = binning.begin(); ib != binning.end(); ++ib)
	for (const auto &ib : _binning)
	{
		//returns only non-flat bins
		if (ib.second.size() < 2)
			continue;

		std::vector<double> bins(ib.second.size() + 1);

		//bin width
		double bw = std::abs(ib.second[1] - ib.second[0]);

		//fill first bin
		bins[0] = ib.second[0] - bw / 2.;

		//and add step until end
		for (size_t i = 1; i < bins.size(); ++i)
			bins[i] = bins[i-1] + bw;

		histbin[ib.first] = bins;
	}

	return histbin;
}

int ParameterSpace::GetEntries()
{
	int npoints = 1;
	for (const auto &ib : _binning)
		npoints *= ib.second.size();

	return npoints;
}

double ParameterSpace::GetPenalty(int n)
{
	std::map<std::string, double> parms = GetEntry(n);

	double x2 = 0;
	for (const auto &ip : _penals) {
		if (ip.second.size() < 2)
			continue;
		x2 += pow((parms[ip.first] - ip.second[0]) / ip.second[1], 2);
	}

	return x2;
}

void ParameterSpace::GetEntry(int n, double &M12, double &M23,
			      double &S12, double &S13, double &S23, double &dCP)
{
	std::map<std::string, double> parms = GetEntry(n);
	M12 = parms["M12"];
	M23 = parms["M23"];
	S12 = parms["S12"];
	S13 = parms["S13"];
	S23 = parms["S23"];
	dCP = parms["CP"];
}

std::map<std::string, double> ParameterSpace::GetEntry(int n)
{
	std::map<std::string, double> vars;

	int a = 0, q = 1;
	Binning::reverse_iterator ir;
	for (ir = _binning.rbegin(); ir != _binning.rend(); ++ir)
	{
		n = (n - a) / q;

		q = ir->second.size();
		a = n % q;

		vars[ir->first] = ir->second.at(a);
	}

	return vars;
}

void ParameterSpace::GetNominal(double &M12, double &M23,
			      double &S12, double &S13, double &S23, double &dCP)
{
	std::map<std::string, double> noms = GetNominal();
	M12 = noms["M12"];
	M23 = noms["M23"];
	S12 = noms["S12"];
	S13 = noms["S13"];
	S23 = noms["S23"];
	dCP = noms["CP"];
}

std::map<std::string, double> ParameterSpace::GetNominal()
{
	std::map<std::string, double> vars;
	for (const auto &ir : _binning)
		vars[ir.first] = ir.second[_nominal[ir.first]];

	return vars;
}

int ParameterSpace::GetNominalEntry()
{
	int n = 0, q = 1;
	Binning::reverse_iterator ir;
	for (ir = _binning.rbegin(); ir != _binning.rend(); ++ir) {
		//std::cout << ir->first << "\t" << ir->second.size() << "\t" << _nominal[ir->first] << std::endl;
		n += _nominal[ir->first] * q;
		q *= ir->second.size();
	}

	return n;
}

std::vector<int> ParameterSpace::GetScanEntries(const std::vector<std::string> &p)
{
	for (const std::string &ip : p)
		if (!_nominal.count(ip)) {
			std::string knowns;
			for (const auto &in : _nominal)
				knowns += " \"" + in.first + "\"";
			throw std::invalid_argument("ParameterSpace:: unknown parameter \""
					+ ip + "\"\n\trecognised parameters are" + knowns);
		}

	std::set<std::string> parms(p.begin(), p.end());

	std::map<std::string, size_t> index;
	std::map<std::string, size_t> maxes;
	for (const auto &in : _nominal) {
		if (parms.find(in.first) == parms.end())
			index[in.first] = in.second;
		else {
			index[in.first] = 0;
			maxes[in.first] = _binning[in.first].size();
		}
	}

	std::vector<int> entries;
	size_t reset = 0;
	do {
		int n = 0, q = 1;
		Binning::reverse_iterator ir;
		for (ir = _binning.rbegin(); ir != _binning.rend(); ++ir) {
			n += index[ir->first] * q;
			q *= ir->second.size();
		}
		entries.push_back(n);

		reset = 0;
		auto ii = index.rbegin();
		for ( ; ii != index.rend(); ++ii) {
			if (parms.find(ii->first) == parms.end())
				continue;

			// index can be incremented
			if (ii->second < _binning[ii->first].size() - 1) {
				++(ii->second);
				break;
			}
			// reset it to zero and move to next
			else {
				ii->second = 0;
				++reset;
			}
		}

	} while (reset < parms.size());

	return entries;
}
