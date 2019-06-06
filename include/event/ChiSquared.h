/* ChiSquared
 * compute the likelihood ratio chi-squared using the pull approach to include systematics
 */



#ifndef ChiSquared_H
#define ChiSquared_H

class ChiSquared
{
	public:
		ChiSquared();

	private:

		//Q matrix is Fij * sigma_j or 
		Eigen::MatrixXd sysQ;
		//numE vector is experiment, numT is expected in theory
		//pull is vecor of systematics pulls that is minimised
		Eigen::VectorXd evtO, evtE, pull, gamma;

};

#endif
