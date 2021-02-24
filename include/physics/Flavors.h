#ifndef FLAVORS_H
#define FLAVORS_H

#include <string>

struct Nu
{
	enum Flavor {
		E_ = 0,
		electron = E_,
		M_ = 1,
		muon = M_,
		T_ = 2,
		tau = T_,
		Eb = 3,
		antielectron = Eb,
		Mb = 4,
		antimuon = Mb,
		Tb = 5,
		antitau = Tb,

		_undefined = -1
	};

	static Flavor fromString(const std::string &flv) {
		if (flv == "nuEb")
			return Flavor::Eb;
		else if (flv == "nuMb")
			return Flavor::Mb;
		else if (flv == "nuTb")
			return Flavor::Tb;
		else if (flv == "nuE0")
			return Flavor::E_;
		else if (flv == "nuM0")
			return Flavor::M_;
		else if (flv == "nuT0")
			return Flavor::T_;
		else
			return _undefined;
	}
	
	static Flavor fromPDG(int pdg) {
		switch (pdg) {
			case -12:
				return Flavor::Eb;
			case -14:
				return Flavor::Mb;
			case -16:
				return Flavor::Tb;
			case 12:
				return Flavor::E_;
			case 14:
				return Flavor::M_;
			case 16:
				return Flavor::T_;
			default:
				return _undefined;
		}
	}

	static std::string toString(int pdg) {
		return toString(fromPDG(pdg));
	}

	static std::string toString(Flavor flv) {
		switch (flv) {
			case Flavor::Eb:
				return "nuEB";
			case Flavor::Mb:
				return "nuMB";
			case Flavor::Tb:
				return "nuTB";
			case Flavor::E_:
				return "nuE0";
			case Flavor::M_:
				return "nuM0";
			case Flavor::T_:
				return "nuT0";
			default:
				return "";
		}
	}
};

#endif
