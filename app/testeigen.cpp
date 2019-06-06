#include <iostream>

#include "Eigen/Dense"

int main()
{
	Eigen::RowVectorXd uu(5);
	uu << 1, 2, 3, 4, 5;

	Eigen::MatrixXd dd = uu.adjoint() * uu;

	std::cout << uu.adjoint() * uu << std::endl;
	std::cout << uu.adjoint() * uu - Eigen::MatrixXd(uu.asDiagonal()) << std::endl;

	return 0;
}
