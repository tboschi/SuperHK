#include "ChiSquared.h"

/* constructor of class
 * requires event stackers of observed events and expected events and file of systematics
 */
ChiSquared::ChiSquared(EventStacker *&obsr, EventStacker *&expt, std::string systFile) :
{
	//load Fij tables and sigma_j from root file
	
	//create matrix Q using Fij and sigma_j
	sysQ = LoadSystematics(systFile);
	//create vector of observed bins
	evtO = StackToVector(obsr);
	//create vector of expected bins
	evtE = StackToVector(expt);
}

Eigen::MatrixXd ChiSquared::LoadSystematics(std::string systFile)
{
	TFile *syst = new TFile(systFile.c_str(), "OPEN");
	return mat;
}

Eigen::VectorXd ChiSquared::StackToVector(EventStacker *& evtStack)
{
	//create a vector with one entry per bin
	std::vector<double> vBin = evtStack->Vectorise();
	Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(vBin.data(), vBin.size());

	return vec;
}

/* to find the minimum of the pulls a system Ax = b must be solved,
 * where x are the pulls.
 * the system has dimension k = number of systematic errors
 * matrix A is returned by the function _A()
 * and the output b is returned by _b()
 */
Eigen::MatrixXd ChiSquared::_A()
{
	return Eigen::MatrixXd::Identity(nSyst) + systQ.transpose() * evtO.asDiagonal() * systQ;
}

Eigen::VectorXd ChiSquared::_b()
{
	//gamma is a vector w/ length n = nBins from Q * pulls
	Eigen::VectorXd g = gamma();
	//if pulls are zero there is no point on computing gamma for b
	if (!g.isZero())
	{
		//add corrections to gamma as in g + 1/(1+g)
		for (unsigned int n = 0; n < gamma.size(); ++n)
			g(n) += 1.0/(1.0 + g(n));
		return systQ.transpose() * (gamma.asDiagonal() * evtO - evtE);
	}
	else
		return systQ.transpose() * (evtO - evtE);
}

double ChiSquared::ComputeX2()
{
	//set initial pulls to zero
	//at each iterations pulls (and gammas) should move away from zero to a stable solution
	pull.setZero();

	double x2 = X2(), x2prev = x2 + 2*conv;
	while ( std::abs(x2 - x2prev) > conv )
	{
		//save previous x2;
		x2prev = x2;

		//pull is the solution of _A * x = _b
		pull = _A().llt().solve(_b());
		//compute chisquared
		x2 = X2();

	}	//while convergence is reached

	return x2;
}

Eigen::VectorXd gamma()
{
	//if pulls is zero don't multiply
	if (!pulls.isZero())
		return systQ * pulls;
	else
		return pulls;
}

double ChiSquared::X2()
{
	return ObservedX2() + SystematicX2();
}

double ChiSquared::ObservedX2()
{
	Eigen::VectorXd gammaVec = gamma();

	double x2 = 0.0;
	for (unsigned int n = 0; i n < nBins; ++n)
	{
		double e = evtE(n);
		double o = evtO(n);
		double g = gammaVec(n);

		x2 += 2 * ( e * (1 + g) - o(n) * (1 + log(1 + g) + log(e / o) ) )
	}

	return x2;
}

double ChiSquared::SystematicX2()
{
	return pull.squaredNorm();
}
