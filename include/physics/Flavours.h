#ifndef Flavours_H
#define Flavours_H

#include <string>

struct Nu
{
	enum Flavour {
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

	static Flavour fromString(const std::string &flv) {
		if (flv == "nuEb")
			return Flavour::Eb;
		else if (flv == "nuMb")
			return Flavour::Mb;
		else if (flv == "nuTb")
			return Flavour::Tb;
		else if (flv == "nuE0")
			return Flavour::E_;
		else if (flv == "nuM0")
			return Flavour::M_;
		else if (flv == "nuT0")
			return Flavour::T_;
		else
			return _undefined;
	}
	
	static Flavour fromPDG(int pdg) {
		switch (pdg) {
			case -12:
				return Flavour::Eb;
			case -14:
				return Flavour::Mb;
			case -16:
				return Flavour::Tb;
			case 12:
				return Flavour::E_;
			case 14:
				return Flavour::M_;
			case 16:
				return Flavour::T_;
			default:
				return _undefined;
		}
	}

	static std::string toString(int pdg) {
		return toString(fromPDG(pdg));
	}

	static std::string toString(Flavour flv) {
		switch (flv) {
			case Flavour::Eb:
				return "nuEB";
			case Flavour::Mb:
				return "nuMB";
			case Flavour::Tb:
				return "nuTB";
			case Flavour::E_:
				return "nuE0";
			case Flavour::M_:
				return "nuM0";
			case Flavour::T_:
				return "nuT0";
			default:
				return "";
		}
	}
};

#endif
