/*
 * Tools
 * Const: standard model constants
 *
 * Author: Tommaso Boschi
 */

#ifndef CONST_H
#define CONST_H

#include <cmath>

//Constants
namespace Const
{
	static const double m  = 1.0;
	static const double km = 1.0e3  * m;
	static const double cm = 1.0e-2 * m;
	static const double mm = 1.0e-3 * m;

	static const double s  = 1.0;

	static const double GeV = 1.0;
	static const double TeV = 1.0e3  * GeV;
	static const double MeV = 1.0e-3 * GeV;
	static const double keV = 1.0e-6 * MeV;
	static const double eV  = 1.0e-9 * keV;

	static const double fC = 299792458 * m / s;		//m/s
	static const double fhBar = 6.5821189916e-25 * GeV * s;	//GeV s, from PDG
	static const double fhBarC = fC * fhBar * GeV * m;
	static const double fAem = 1.0/137.035999074;		// from PDG
	static const double fNa = 6.02214085774e23;	//mol-1

	//Conversion
	static const double fM2GeV = 5.06e15;		//1GeV in 1/m
	static const double fS2GeV = 1.52e24;		//1GeV in 1/s
	static const double fGeV2cm = 389.4e-30;	//1GeV-2 in cm2
	static const double fGeV2ub = 0.3894e3;		//1GeV-2 in ub
	static const double fPi = 3.1415926536;		//pi
	static const double fPi2 = fPi*fPi;		//pi
	static const double fPi3 = fPi2*fPi;		//pi
	static const double fPi4 = fPi3*fPi;		//pi
	static const double fPi5 = fPi4*fPi;		//pi
	static const double fDeg = 180.0/fPi;		//Rad to Deg

	//SM constant - PDG 2016
	static const double fGF = 1.16637876e-5 / GeV / GeV;	//Fermi constant in natural units
	static const double fGF2 = fGF*fGF;			//From PDG
	static const double fSin2W = 0.23129;		//Sin weinberg squared - MSbar scheme
	static const double fL2E = 2.534;		//L km over E GeV conversion

	//CKM entries
	static const double fU_ud = 0.97417;
	static const double fU_us = 0.2248;
	static const double fU_ub = 0.0409;
	static const double fU_cd = 0.220;
	static const double fU_cs = 0.995;
	static const double fU_cb = 0.0405;
	static const double fU_td = 0.0082;
	static const double fU_ts = 0.0400;
	static const double fU_tb = 1.009;

	//PMNS entries
	static const double fU_e1 = 0.81;
	static const double fU_e2 = 0.54;
	static const double fU_e3 = -0.15;
	static const double fU_m1 = -0.35;
	static const double fU_m2 = 0.70;
	static const double fU_m3 = 0.62;
	static const double fU_t1 = 0.44;
	static const double fU_t2 = -0.45;
	static const double fU_t3 = 0.77;

	//Masses in GeV - PDG 2017
	static const double fMQuarkU   = 2.2e-3 * GeV;
	static const double fMQuarkD   = 4.7e-3 * GeV;
	static const double fMQuarkS   = 96e-3 * GeV;
	static const double fMQuarkC   = 1.28 * GeV;
	static const double fMQuarkB   = 4.18 * GeV;
	static const double fMQuarkT   = 173.1 * GeV;
	static const double fMElectron = 0.510999e-3 * GeV;
	static const double fMMuon     = 105.6583e-3 * GeV;
	static const double fMTau      = 1776.86e-3 * GeV;
	static const double fMPion     = 139.57061e-3 * GeV;
	static const double fMPion0    = 134.9770e-3 * GeV;
	static const double fMKaon     = 493.677e-3 * GeV;
	static const double fMKaon0    = 497.611e-3 * GeV;
	static const double fMEta      = 547.862e-3 * GeV;
	static const double fMEtai     = 957.78e-3 * GeV;
	static const double fMRho      = 775.11e-3 * GeV;
	static const double fMRho0     = 775.26e-3 * GeV;
	static const double fMOmega    = 782.65e-3 * GeV;
	static const double fMKaonx    = 891.76e-3 * GeV;
	static const double fMKaon0x   = 682e-3 * GeV;
	static const double fMPhi      = 1019.460e-3 * GeV;
	static const double fMCharm    = 1869.56e-3 * GeV;
	static const double fMDs       = 1968.28e-3 * GeV;
	static const double fMProton   = 938.272081e-3 * GeV;
	static const double fMNeutron  = 939.565143e-3 * GeV;
	static const double fMW        = 80.385 * GeV;
	static const double fMZ        = 91.1876 * GeV;

	
	//decay constant - 0901.3789
	/*
	static const double fDPion   = 0.1327;
	static const double fDPion0  = 0.1300;
	static const double fDKaon   = 0.1598;
	static const double fDRho    = 0.2200;
	static const double fDRho0   = 0.2200;
	static const double fDEta    = 0.1647;
	static const double fDEtai   = 0.1529;
	static const double fDOmega  = 0.1950;
	static const double fDKaonx  = 0.2170;
	static const double fDKaon0x = 0.2170;	
	static const double fDPhi    = 0.2290;
	static const double fDCharm  = 0.2226;
	*/

	//decay constant - 1805.08567
	static const double fDPion   = 0.1302;		//GeV
	//static const double fDPion0  = 0.1300;	//GeV
	static const double fDKaon   = 0.1556;		//GeV
	static const double fDRho    = 0.1620;		//GeV²
	static const double fVRho    = 1 - 2 * fSin2W;
	//static const double fDRho0   = 0.2200;	//GeV²
	static const double fDEta    = 0.0817;		//GeV
	static const double fDEtai   = 0.0947;		//GeV
	static const double fDOmega  = 0.1530;		//GeV²
	static const double fVOmega  = 4 * fSin2W / 3.0;
	static const double fDKaonx  = 0.1933;	//maybe not possible
	//static const double fDKaon0x = 0.2170;	
	static const double fDPhi    = 0.2340;		//GeV²
	static const double fVPhi    = 4 * fSin2W / 3.0 - 1;
	static const double fDCharm  = 0.2120;		//GeV

	static const double fKaPi = 0.9700;	//f(0)
	static const double fK0L_ = 0.0267;	//Linear dependence of f+ in K0m3 (PDG)
	static const double fK0L0 = 0.0117;	//Linear dependence of f0 in K0m3 (PDG)
	static const double fKCL_ = 0.0277;	//Linear dependence of f+ in K+m3 (PDG)
	static const double fKCL0 = 0.0183;	//Linear dependence of f0 in K+m3 (PDG)

	static const double fMagMuN = -1.9130427345;	//neutron magnetic moment (in nuclear magneton units)
	static const double fMagMuP = 2.79284735128;	//proton magnetic moment (in nuclear magneton units);
	static const double fMA = 0.990;		//GeV, axial mass, from GENIE
	//static const double fMA = 1.032;		//GeV, axial mass, from PRD35, 785 (1987) [Arhens]
	//static const double fMA = 1.270;		//GeV, axial mass, from PRD92, 113011 (2015) [lattice]
	//static const double fMA = 1.026;		//GeV, axial mass, from Giunti-Kim
	static const double fMV = 0.840;		//GeV, vectorial mass, from GENIE
	static const double fGA0 = 1.2671;		//Axial form factor at Q2 = 0, from GENIE

	static const double fWwidth = 2.085;		//W decay width, in GeV
	static const double fZwidth = 2.4952;		//W decay width, in GeV
}

#endif
