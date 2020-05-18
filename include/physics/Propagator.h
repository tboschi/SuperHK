/* namespace continaing some functions for Oscillator class
 * not sure if makes sense to keep here
 * rather than moving into Oscillator class
 */

#ifndef Propagator_H
#define Propagator_H

namespace Osc
{
	static Eigen::VectorXd MatterStates(const Eigen::MatrixXcd &pmnsM,	//PMNS matrix
					    const Eigen::VectorXd &dms,		//masses in vacuum
					    const double &ff)	       		 //density factor
	{
		double ms1   =  dms(0);		//mass squared
		double dms12 = -dms(1);		//delta m squared
		double dms13 = -dms(2);		//delta m squared

		double Ue1s  = std::norm(pmnsM(0, 1));	//Ue1 squared
		double Ue2s  = std::norm(pmnsM(0, 2));	//Ue1 squared

		double alpha = ff + dms12 + dms13; 

		double beta = dms12 * dms13 + ff * (dms12 * (1 - Ue2s) + dms13 * (1 - Ue1s) );

		double gamma = ff * dms12 * dms13 * Ue1s;
		
		double argo = (2 * pow(alpha, 3) - 9 * alpha * beta + 27 * gamma) /
			      pow( pow(alpha, 2) - 3 * beta, 1.5) / 2.0;

		if (std::abs(argo) > 1.0)
			argo /= std::abs(argo);

		double r0 = acos(arg) / 3.0;
		double r1 = root0 - 2.0/3.0 * Const::fPi;
		double r2 = root0 + 2.0/3.0 * Const::fPi;

		Eigen::VectorXd vMat(3);
	  
		vMat(0) = - 2.0/3.0 * sqrt(alpha * alpha - 3 * beta) * cos(r0) + ms1 - alpha/3.0;
		vMat(1) = - 2.0/3.0 * sqrt(alpha * alpha - 3 * beta) * cos(r1) + ms1 - alpha/3.0;
		vMat(2) = - 2.0/3.0 * sqrt(alpha * alpha - 3 * beta) * cos(r2) + ms1 - alpha/3.0;

		return vMat;
	}

	//This is equivalent to get_product in mosc.cc
	//works with any neutrino dimension
	static Eigen::Matrixcd TransitionMatrix(const Eigen::MatrixXcd &dmMatVac, //masses
						const Eigen::MatrixXcd &dmMatMat, //masses
						const Eigen::RowVectorXcd &Ue,	//PMNS elec row
						const double &ff,		//density factor
						const double &l2e)		//L over E factor

	{
		unsigned int nDim = Ue.size();	//number of columns

		Eigen::MatrixXcd product(nDim);
		Eigen::MatrixXcd result(nDim);

		for (unsigned int k = 0; k < nDim; ++k)
		{
			product.setIdentity();		//set product to I(nDim)
			for (unsigned int j = 0; j < nDim; ++j)
			{
				if (j == k)		//product for j != k
					continue;

				//this is (2EH-M / dM²)_j
				product *= (-ff * Ue.conjugate() * Ue - 
						dmMatVac.row(j).asDiagonal()) / 
					   dmMatMat(k, j);
			}

			double arg = - dmMatVac(k, 0) * l2e;
			product *= std::complex<double>(cos(arg), sin(arg));

			result += product;
		}

		return result;
	}
}

#endif
