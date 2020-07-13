#ifndef CONST_H
#define CONST_H

//Constants
namespace Const
{
	static const double C = 299792458;		//m/s
	static const double hBar = 6.5821189916e-25;	//GeV s, from PDG
	static const double hBarC = C * hBar;		// GeV * m
	static const double Aem = 1.0/137.035999074;	// from PDG
	static const double Na = 6.02214085774e23;	//mol-1

	static const double GF = 1.16637876e-5;		//Fermi constant in natural units
	static const double GF2 = GF*GF;		//From PDG

	static const double L2E = 2.534;		//L km over E GeV conversion

	static const double pi = 3.1415926536;		//pi
	static const double pi2 = pi*pi;		//pi
	static const double pi3 = pi2*pi;		//pi
	static const double pi4 = pi3*pi;		//pi
	static const double pi5 = pi4*pi;		//pi
	static const double Deg = 180.0/pi;		//Rad to Deg
}

#endif
