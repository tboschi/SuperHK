#include "physics/Oscillator.h"

Oscillator::Oscillator(const std::vector<double> &lengths,
		       const std::vector<double> &densities,
		       int numDimension, bool lut) :
	nDim(numDimension),
	kLUT(lut)
{
	for (int i = 0; i < lengths.size(); ++i)
		lens_dens.push_back(std::make_pair(lengths[i], densities[i]));

	transM = Eigen::MatrixXcd::Identity(nDim, nDim);
}

Oscillator::Oscillator(CardDealer *cd, int numDimension) :
	nDim(numDimension)
{
	std::string densityFile;
	if (cd->Get("density_profile", densityFile))
		SetDensity(densityFile);
	else {	//default  vacuum oscillation
		lens_dens.clear();
		lens_dens.push_back(std::make_pair(295., 0));
	}

	int lut;
	if (!cd->Get("LUT", lut))	//look up table stores matrices
		kLUT = false;
	else
		kLUT = lut;

	transM = Eigen::MatrixXcd::Identity(nDim, nDim);
}

Oscillator::Oscillator(const std::string &densityFile, int numDimension) :
	nDim(numDimension)
{
	SetDensity(densityFile);

	transM = Eigen::MatrixXcd::Identity(nDim, nDim);
}

void Oscillator::SetDensity(const std::string &densityFile)
{
	lens_dens.clear();

	std::ifstream inf(densityFile);
	std::string line;
	double ll, dd;
	//std::cout << "Reading " << densityFile << std::endl;
	while(std::getline(inf, line))
	{
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));

		if (line.empty())
			continue;

		std::stringstream ssl(line);
		if (ssl >> ll >> dd)
			lens_dens.push_back(std::make_pair(ll, dd));
	}
}

//return oscillation probability from falvour in to flavour out at energy given
//the lengths and densities are fixed beforehand
double Oscillator::Probability(Nu in, Nu out, double energy)
{
	if (in * out < 0)
	{
		std::cerr << "OScillator::Probability can't have oscillation from neutrino to antineutrino\n"
			  << "			      and vice versa!" << std::endl;
		return -1.0;
	}

	kAntineutrino = (in < 0);				//flag if antineutrino
	double ten = energy * (kAntineutrino ? -1.0 : 1.0);	//test energy, negative for antineutr

	int iFlav = std::abs(in)  - 1;
	int oFlav = std::abs(out) - 1;

	//transposition done here!! should save some little copmutation time
	//int nEntr = oFlav + iFlav*nDim;
	Eigen::MatrixXd tM = TransitionMatrix(energy).cwiseAbs2();
	//std::vector<double> tV(tM.data(), tM.data() + tM.size());
	//std::cout << "Matrix " << tM(oFlav, iFlav) << "\tvector " << tV.at(nEntr) << std::endl;
	//std::cout << "Matrix\n" << tM << std::endl << std::endl;
	//return tV.at(nEntr);
	return tM(oFlav, iFlav);

	auto ilut = FindEnergy(ten);

	if (ilut != mLUT.end())	//value already computed and save
		return (ilut->second)(oFlav, iFlav);
	else			//save computation in LUT
	{
		mLUT[ten] = TransitionMatrix(energy).cwiseAbs2();

		return mLUT[ten](oFlav, iFlav);
	}
}


void Oscillator::Oscillate(Nu in, Nu out, TH1D* h)
{
	for (int i = 1; i < h->GetNbinsX()+1; ++i)
	{
		double en = h->GetBinCenter(i);
		double cn = h->GetBinContent(i);

		//if (std::abs(cn) < 1e-9)
		//	continue;

		h->SetBinContent(i, cn * Probability(in, out, en));
		//std::cout << "bin " << en << " : " << cn << " -> " << h->GetBinContent(i) << std::endl;
	}
}

void Oscillator::Reset()
{
	mLUT.clear();
}

std::map<double, Eigen::MatrixXd>::iterator Oscillator::FindEnergy(double energy)
{
	auto ilut = mLUT.lower_bound(energy);

	if (ilut == mLUT.end())
		return ilut;

	if (std::abs(ilut->first - energy) < 1e-9)
		return ilut;

	if (ilut != mLUT.begin() && std::abs(std::prev(ilut)->first - energy) < 1e-9)
		return std::prev(ilut);

	return mLUT.end();
}

double Oscillator::Length() {
	return std::accumulate(lens_dens.begin(), lens_dens.end(), 0,
			[&](double sum, std::pair<double, double> ild)
			{ return sum + ild.first; });
}

double Oscillator::Density() {
	return std::accumulate(lens_dens.begin(), lens_dens.end(), 0,
			[&](double sum, std::pair<double, double> ild)
			{ return sum + ild.first * ild.second; }) / Length();
}

/////////////**********************ENTERING THE WORLD OF PHYSICS AND BLACK MAGIC HERE
//
//
//
//
//This is equivalent to propagate(int) in BargerPropagator.cc
Eigen::MatrixXcd Oscillator::TransitionMatrix(double energy)
{
	transM.setIdentity();

	if (kAntineutrino)
		pmnsM = _pmnsM.conjugate();
	else
		pmnsM = _pmnsM;

	int sign = 2*kAntineutrino - 1;			//+1 if antineutrino, -1 if neutrino
	for (const auto &ld : lens_dens)
	{
		//std::cout << "Enu " << energy << "\trho " << vDensity.at(i)/2.0 << std::endl;
		//std::cout << "sqrtFG " << sqrt(8) * fG << "\tsign " << sign << std::endl;
		double ff  = sign * sqrt(8) * fG * energy * ld.second/2.0;
		double l2e = Const::fL2E * ld.first/energy;

		//std::cout << "factor " << ff << std::endl;

		transM *= TransitionMatrix(ff, l2e);
	}

	//if (kAntineutrino)
	//	transM = pmnsM.conjugate() * transM * pmnsM.transpose();
	//else
	//	transM = pmnsM * transM * pmnsM.adjoint();
	//std::cout << "\nmix\n" << transM.cwiseAbs2() << std::endl;

	transM = pmnsM * transM * pmnsM.adjoint();
	return transM;
}

//This is equivalent to getA in mosc.cc
Eigen::MatrixXcd Oscillator::TransitionMatrix(double ff, double l2e)
{
	Eigen::MatrixXd dmMatVac(nDim, nDim), dmMatMat(nDim, nDim);
	Eigen::MatrixXcd product(nDim, nDim),
			 result = Eigen::MatrixXcd::Zero(nDim, nDim);
	Eigen::RowVectorXcd Ue = pmnsM.row(0);

	//load matter matrices
	MatterMatrices(dmMatVac, dmMatMat, ff);

	//std::cout << "\nMat\n" << dmMatMat << "\nVac\n" << dmMatVac << std::endl;
	for (int k = 0; k < nDim; ++k)
	{
		product.setIdentity();		//set product to I(nDim)
		for (int j = 0; j < nDim; ++j)
		{
			if (j == k)		//product for j != k
				continue;

			//this is (2EH-M / dM²)_j
			product *= (-ff * Ue.adjoint() * Ue - 
				    Eigen::MatrixXd(dmMatVac.row(j).asDiagonal()) ) / 
				    dmMatMat(k, j);
		}

		double arg = - dmMatVac(k, 0) * l2e;
		product *= std::complex<double>(std::cos(arg), std::sin(arg));

		result += product;
	}

	return result;
}

//This is equialent to getM in mosc.cc
// The strategy to sort out the three roots is to compute the vacuum
// mass the same way as the "matter" masses are computed then these are sorted
// according to the vaccum solutions
void Oscillator::MatterMatrices(Eigen::MatrixXd &dmMatVac,	//output - mass diff matter-vacuum
				Eigen::MatrixXd &dmMatMat,	//output - mass diff matter
				const double &ff)		//density factor
{
	Eigen::VectorXd vVac = MatterStates(0.0);	//vacuum solutions
	Eigen::VectorXd vMat = MatterStates(ff);	//matter solutions
  
	////std::cout << "vVac\n" << vVac << std::endl << std::endl;
	//std::cout << "vMat\n" << vMat << std::endl << std::endl;
	//sorting according to which matter solution is closest to the respective vacuum sol.
	//
	for (int i = 0; i < nDim-1; ++i)
	{
		int k = i;
		double val = fabs(dms(i, 0) - vVac(i));
		for (int j = i+1; j < nDim; ++j)
		{
			if (val > fabs(dms(i, 0) - vVac(j)))
			{
				k = j;
				val = fabs(dms(i, 0) - vVac(j));
			}
		}
		//swap -- not sure if std::swap would work here
		//val = vMat(i);
		//vMat(i) = vMat(k);
		//vMat(k) = val;
		std::swap(vMat(i), vMat(k));
	}

	Eigen::MatrixXd mVac(nDim, nDim), mMat(nDim, nDim);
	//matrix made ouf of vector Mat
	mMat << vMat, vMat, vMat;	/*	M1²	M1²	M1²
						M2²	M2²	M2²
						M3²	M3²	M3² */

	//matrix made ouf of vector vacuum mass differences
	mVac << dms, dms, dms;		/*	m1²-m1²,m1²-m1²,m1²-m1²
						m2²-m1²,m2²-m1²,m2²-m1²
						m3²-m1²,m3²-m1²,m3²-m1² */

	dmMatMat = (mMat - mMat.transpose());	/*	M1²-M1²,M1²-M2²,M1²-M3²
							M2²-M1²,M2²-M2²,M2²-M3²
							M3²-M1²,M3²-M2²,M3²-M3² */

	dmMatVac = (mMat - mVac.transpose());	//mixture of the above

	//std::cout << "dmMatMat\n" << dmMatMat << std::endl << std::endl;
	//std::cout << "dmMatVac\n" << dmMatVac << std::endl << std::endl;

}

// compute the matter mass squared vector, solution of Eq. 21/22 if PRD 22.11 (1980)
// the input vacuum mass squared vector is defined as (m1², m2²-m1², m3²-m1²)
//works only with 3 neutrino states
Eigen::VectorXd Oscillator::MatterStates(const double &ff)	//density factor	
{
	//pmns and massSquare are global
	double ms1   =  dms(0);		//mass squared
	double dms12 = -dms(1);		//delta m squared
	double dms13 = -dms(2);		//delta m squared

	double Ue1s  = std::norm(pmnsM(0, 0));	//Ue1 squared
	double Ue2s  = std::norm(pmnsM(0, 1));	//Ue2 squared
	double Ue3s  = std::norm(pmnsM(0, 2));	//Ue3 squared

	double alpha = ff + dms12 + dms13; 

	double beta  = dms12 * dms13 + ff * (dms12 * (1 - Ue2s) + dms13 * (1 - Ue3s) );

	double gamma = ff * dms12 * dms13 * Ue1s;
	
	double argo  = (alpha * (2 * pow(alpha, 2) - 9 * beta) + 27 * gamma) /
		       pow( pow(alpha, 2) - 3 * beta, 1.5) / 2.0;

	if (std::abs(argo) > 1.0)
		argo /= std::abs(argo);

	double r0 = std::acos(argo) / 3.0;
	double r1 = r0 - 2.0/3.0 * Const::fPi;
	double r2 = r0 + 2.0/3.0 * Const::fPi;

	//std::cout << "argo " << argo << std::endl;

	Eigen::VectorXd vMat(3);
  
	vMat(0) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta) * std::cos(r0) + ms1 - alpha/3.0;
	vMat(1) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta) * std::cos(r1) + ms1 - alpha/3.0;
	vMat(2) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta) * std::cos(r2) + ms1 - alpha/3.0;

	return vMat;
}

//set pmns and masses
void Oscillator::SetMasses_NH(double dms21, double dms23)
{
	dms = Eigen::VectorXd::Zero(nDim);
	dms << 0, dms21, dms21+dms23;
}

void Oscillator::SetMasses_IH(double dms21, double dms23)
{
	SetMasses_NH(dms21, -dms23-dms21);
}

void Oscillator::SetMasses_abs(double ms2, double ms3)
{
	dms = Eigen::VectorXd::Zero(nDim);
	dms << 0, ms2, ms3;
}

void Oscillator::SetPMNS_sin(double s12, double s13, double s23, double cp) 
{
	_pmnsM = Eigen::MatrixXd::Zero(nDim, nDim);

	double c12 = sqrt(1 - s12*s12);
	double c13 = sqrt(1 - s13*s13);
	double c23 = sqrt(1 - s23*s23);

	//exp(-i cp)
	std::complex<double> dcp(std::cos(cp), -std::sin(cp));

	Eigen::Matrix3cd U1, U2, U3;

	U1 <<	1.0, 	0.0, 	0.0,
	   	0.0, 	c23,	s23,
		0.0, 	-s23, 	c23;

	U2 <<	c13,		0.0, 	s13*dcp,
	   	0.0, 		1.0,	0.0,
		-s13/dcp, 	0.0,	c13;

	U3 <<	c12, 	s12,	0.0,
	   	-s12, 	c12,	0.0,
		0.0,	0.0,	1.0;

	_pmnsM = U1 * U2 * U3;
}

//passing sin squared
void Oscillator::SetPMNS_sin2(double s12, double s13, double s23, double cp) 
{
	SetPMNS_sin(sqrt(s12), sqrt(s13), sqrt(s23), cp);
}

//passing angles in rad
void Oscillator::SetPMNS_angles(double t12, double t13, double t23, double cp) 
{
	SetPMNS_sin(std::sin(t12), std::sin(t13), std::sin(t23), cp);
}

Eigen::MatrixXcd Oscillator::PMNS()
{
	return _pmnsM;
}

Eigen::VectorXd Oscillator::Masses()
{
	return dms;
}
