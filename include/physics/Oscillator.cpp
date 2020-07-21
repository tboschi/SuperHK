#include "physics/Oscillator.h"

Oscillator::Oscillator(const std::vector<double> &lengths,
		       const std::vector<double> &densities,
		       int neutrino, bool lut, double threshold) :
	fG(Const::GF * Const::Na * pow(Const::hBarC * 1e8, 3)),
	_dim(neutrino),
	_thr(threshold),
	kLUT(lut)
{
	for (int i = 0; i < lengths.size(); ++i)
		lens_dens.push_back(std::make_pair(lengths[i], densities[i]));

	//trans = Eigen::MatrixXcd::Identity(_dim, _dim);
}

Oscillator::Oscillator(const std::string &densityFile, 
		       int neutrino, bool lut, double threshold) :
	fG(Const::GF * Const::Na * pow(Const::hBarC * 1e8, 3)),
	_dim(neutrino),
	_thr(threshold),
	kLUT(lut)
{
	SetMatterProfile(densityFile);
}



Oscillator::Oscillator(CardDealer *cd) :
	fG(Const::GF * Const::Na * pow(Const::hBarC * 1e8, 3))
{
	std::string densityFile;
	if (cd->Get("density_profile", densityFile))
		SetMatterProfile(densityFile);
	else {	//default  vacuum oscillation
		lens_dens.clear();
		lens_dens.push_back(std::make_pair(295., 0));
	}

	if (!cd->Get("neutrinos", _dim))//look up table stores matrices
		_dim = 3;		// default neutrinos

	if (!cd->Get("threshold", _thr))
		_thr = 1e-9;	// default value

	int lut;
	if (!cd->Get("LUT", lut))	//look up table stores matrices
		kLUT = false;
	else
		kLUT = lut;
}

void Oscillator::SetMatterProfile(const std::string &densityFile)
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

void Oscillator::SetMatterProfile(const std::vector<std::pair<double, double> > &l_d)
{
	lens_dens = l_d;
}

//return oscillation probability from falvour in to flavour out at energy given
//the lengths and densities are fixed beforehand
double Oscillator::Probability(Nu in, Nu out, double energy, bool force)
{
	if (std::abs(in - out) >= 3 && !force)
	{
		std::cerr << "WARNING - Oscillator : can't have oscillation "
			  << " between neutrinos and antineutrinos\n"
			  << "        call Oscillator::Probability(..., force = true)\n"
			  << "        to suppress this message"
			  << std::endl;
	}

	if (!kLUT)
		//return TransitionMatrix(energy).cwiseAbs2()(ofl, ifl);
		return TransitionMatrix(energy).cwiseAbs2()(out, in);

	auto ilut = FindEnergy(energy);
	if (ilut != mLUT.end())	// use precomputed matrix
		return (ilut->second)(out, in);
	else	// save new matrix
	{
		mLUT[energy] = TransitionMatrix(energy).cwiseAbs2();
		return mLUT[energy](out, in);
	}
}


Eigen::VectorXd Oscillator::Oscillate(Nu in, Nu out,
				const std::vector<double> &bins)
{
	Eigen::VectorXd vb(bins.size() - 1);
	for (int i = 0; i < vb.size(); ++i)
		vb(i) = Probability(in, out, (bins[i+1] + bins[i]) / 2.);
	return vb;
}

void Oscillator::Reset()
{
	mLUT.clear();
}

std::map<double, Eigen::MatrixXd>::iterator Oscillator::FindEnergy(double energy)
{
	if (!mLUT.size())	// empty, just return end
		return mLUT.end();

	auto ilut = mLUT.lower_bound(energy);

	// there is a lower bound, so check if it is pointing at the right element
	if (ilut != mLUT.end() && std::abs(ilut->first - energy) < _thr)
		return ilut;

	// works even if ilut == mLUT.end(), check previous element which can be good 
	if (ilut != mLUT.begin() && std::abs(std::prev(ilut)->first - energy) < _thr)
		return std::prev(ilut);

	// nothing matches, so return end
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

//*****************************************************************
//	ENTERING THE WORLD OF PHYSICS AND BLACK MAGIC HERE
//*****************************************************************
//
//
//This is equivalent to propagate(int) in BargerPropagator.cc
Eigen::MatrixXcd Oscillator::TransitionMatrix(double energy, double phase)
{
	//trans.setIdentity();
	Eigen::MatrixXcd trans = Eigen::MatrixXcd::Identity(2*_dim, 2*_dim);

	// looping through different layers and densities
	for (auto ld = lens_dens.begin(); ld != lens_dens.end(); ++ld)
	{
		// compute some energy factors...?
		double ff  = -sqrt(8) * fG * energy * ld->second/2.0;
		double l2e = Const::L2E * ld->first/energy;

		// apply phase offset only to penultimate element..why?
		// it is needed for mantle averaging: by default phase = 0
		// apart from mantle avg
		if (ld == lens_dens.end()-2)
			trans *= TransitionMatrix(ff, l2e, phase);
		else
			trans *= TransitionMatrix(ff, l2e);
	}

	return _pmns * trans * _pmns.adjoint();
}

//This is equivalent to getA in mosc.cc
Eigen::MatrixXcd Oscillator::TransitionMatrix(double ff, double l2e, double phase)
{
	Eigen::MatrixXd dmMatVac(2*_dim, _dim),
			dmMatAnt(2*_dim, _dim);

	Eigen::MatrixXcd Ue2 = ff * _pmns.row(_dim).adjoint() * _pmns.row(_dim)
 			     - ff * _pmns.row(0).adjoint() * _pmns.row(0);

	//load matter matrices
	MatterMatrices(dmMatVac, dmMatAnt, ff);

	std::vector<Eigen::MatrixXcd> vmat(_dim,
				      Eigen::MatrixXcd::Identity(2*_dim, 2*_dim));

	Eigen::VectorXd mat_phi = -dmMatVac.row(0) * l2e;
	Eigen::VectorXd ant_phi = -dmMatVac.row(_dim) * l2e;
	mat_phi(_dim-1) += phase;
	ant_phi(_dim-1) += phase;
	for (int k = 0; k < _dim; ++k) {
		vmat[k].topLeftCorner(_dim, _dim) *= std::complex<double>
				(std::cos(mat_phi(k)), std::sin(mat_phi(k)));
		vmat[k].bottomRightCorner(_dim, _dim) *= std::complex<double>
				(std::cos(ant_phi(k)), std::sin(ant_phi(k)));
	}

	for (int j = 0; j < _dim; ++j) {
		//this is (2EH-M / dM²)_j
		Eigen::MatrixXcd eh = Ue2;
		eh.diagonal() -= dmMatVac.col(j);

		for (int k = 0; k < _dim; ++k) {
			if (k == j)
				continue;
			vmat[k] *= eh * (dmMatAnt.col(j) - dmMatAnt.col(k))
					.cwiseInverse().asDiagonal();
		}
	}

	Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(2*_dim, 2*_dim);
	return std::accumulate(vmat.begin(), vmat.end(), result);
}

//This is equialent to getM in mosc.cc
// The strategy to sort out the three roots is to compute the vacuum
// mass the same way as the "matter" masses are computed then these are sorted
// according to the vaccum solutions
void Oscillator::MatterMatrices(Eigen::MatrixXd &dmMatVac, //output - mass diff matter-vacuum
				Eigen::MatrixXd &dmMatAnt, //output - mass diff matter
				double ff)		//density factor
{
	Eigen::VectorXd vVac = MatterStates(0.0);	//vacuum solutions

	Eigen::VectorXd vMat = MatterStates( ff);	//matter solutions
	Eigen::VectorXd vAnt = MatterStates(-ff, _dim);//antimatter solutions

	////std::cout << "vVac\n" << vVac << std::endl << std::endl;
	//std::cout << "vMat\n" << vMat << std::endl << std::endl;
	//sorting according to which matter solution is closest to the respective vacuum sol.
	//
	for (int i = 0; i < _dim-1; ++i) {
		int k = i;
		double val = fabs(dms(i, 0) - vVac(i));
		for (int j = i+1; j < _dim; ++j) {
			if (val > fabs(dms(i, 0) - vVac(j)))
			{
				k = j;
				val = fabs(dms(i, 0) - vVac(j));
			}
		}
		std::swap(vMat(i), vMat(k));
		std::swap(vAnt(i), vAnt(k));
	}

	//matrix made ouf of vector Mat
	dmMatAnt << vMat.transpose(), vMat.transpose(), vMat.transpose(),
		    vAnt.transpose(), vAnt.transpose(), vAnt.transpose();
	//	M1²	M2²	M3²	->  -  1-0 2-0
	//	M1²	M2²	M3²	-> 0-1  -  2-1
	//	M1²	M2²	M3²	-> 0-2 1-2  -
	
	//matrix made ouf of vector vacuum mass differences
	dmMatVac << dms, dms, dms,	//	m1²-m1²,m1²-m1²,m1²-m1²
		    dms, dms, dms;	//	m2²-m1²,m2²-m1²,m2²-m1²
					//	m3²-m1²,m3²-m1²,m3²-m1²

	dmMatVac = dmMatAnt - dmMatVac;
	//Eigen::MatrixXd dmma = dmMatAnt;
	//for (int i = 0; i < _dim; ++i)
	//	dmMatAnt.col(i) -= dmma.col((i + 1) % _dim);

	//std::cout << "newMat\n" << dmMatAnt.topLeftCorner(_dim, _dim) << "\n" << std::endl;
	//std::cout << "newVac\n" << dmMatVac.topLeftCorner(_dim, _dim) << "\n" << std::endl;
}

// compute the matter mass squared vector, solution of Eq. 21/22 if PRD 22.11 (1980)
// the input vacuum mass squared vector is defined as (m1², m2²-m1², m3²-m1²)
//works only with 3 neutrino states
Eigen::VectorXd Oscillator::MatterStates(double ff, int off)	//density factor	
{
	//pmns and massSquare are global
	double ms1   =  dms(0);	//mass squared
	double dms12 = -dms(1);	//delta m squared
	double dms13 = -dms(2);	//delta m squared

	double Ue1s = std::norm(_pmns(off, 0 + off));	//Ue1 squared
	double Ue2s = std::norm(_pmns(off, 1 + off));	//Ue2 squared
	double Ue3s = std::norm(_pmns(off, 2 + off));	//Ue3 squared

	double alpha = ff + dms12 + dms13; 
	double beta = dms12 * dms13 + ff * (dms12 * (1 - Ue2s)
			+ dms13 * (1 - Ue3s) );
	double gamma = ff * dms12 * dms13 * Ue1s;
	
	double argo = (alpha * (2 * pow(alpha, 2) - 9 * beta) + 27 * gamma)
		/ pow( pow(alpha, 2) - 3 * beta, 1.5) / 2.0;

	if (std::abs(argo) > 1.0)
		argo /= std::abs(argo);

	//double r0 = std::acos(argo) / 3.0;
	//double r1 = r0 - 2.0/3.0 * Const::pi;
	//double r2 = r0 + 2.0/3.0 * Const::pi;
	double ang = std::acos(argo) / 3.;

	Eigen::VectorXd vMat(_dim);
	vMat(0) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta)
		* std::cos(ang) + ms1 - alpha/3.0;
	vMat(1) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta)
		* std::cos(ang - 2./3. * Const::pi) + ms1 - alpha/3.0;
	vMat(2) = - 2.0/3.0 * sqrt(pow(alpha, 2) - 3 * beta)
		* std::cos(ang + 2./3. * Const::pi) + ms1 - alpha/3.0;

	return vMat;
}

//set pmns and masses
void Oscillator::SetMasses_NH(double dms21, double dms23)
{
	dms = Eigen::VectorXd::Zero(_dim);
	dms << 0, dms21, dms21+dms23;
}

void Oscillator::SetMasses_IH(double dms21, double dms23)
{
	SetMasses_NH(dms21, -dms23-dms21);
}

void Oscillator::SetMasses_abs(double ms2, double ms3)
{
	dms = Eigen::VectorXd::Zero(_dim);
	dms << 0, ms2, ms3;
}

void Oscillator::SetPMNS_sin(double s12, double s13, double s23, double cp) 
{
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

	_pmns = Eigen::MatrixXcd::Zero(2*_dim, 2*_dim);

	// pmns matrix is 2*ndim X 2*ndim and contains both pmns and pmns*
	_pmns.topLeftCorner(_dim, _dim) = U1 * U2 * U3;
	_pmns.bottomRightCorner(_dim, _dim) =
		_pmns.topLeftCorner(_dim, _dim).conjugate();
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
	return _pmns.topLeftCorner(_dim, _dim);
}

Eigen::VectorXd Oscillator::Masses()
{
	return dms;
}
