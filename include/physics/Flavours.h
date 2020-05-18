#ifndef Flavours_H
#define Flavours_H

enum Nu
{
	//neutrinos
	E_ = 1,
	electron = E_,
	M_ = 2,
	muon = M_,
	T_ = 3,
	tau = T_,
	////////////////
	//antineutrinos
	//they should be all negative
	Eb = -E_,
	antielectron = Eb,
	Mb = -M_,
	antimuon = Mb,
	Tb = -T_,
	antitau = Tb,
};

#endif
