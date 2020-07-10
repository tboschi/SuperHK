#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/CardDealer.h"
#include "physics/Oscillator.h"

#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TKey.h"
#include "TList.h"

int main(int argc, char** argv)
{
	std::vector<TMatrixT<double>*> matrices;
	std::cout << "Reading " << argc-1 << " files" << std::endl;
	for (int f = 1; f < argc; ++f)
	{
		TFile inf(argv[f], "OPEN");
		if (inf.IsZombie())
			continue;
		std::cout << "opening " << argv[f] << std::endl;

		TList *list = inf.GetListOfKeys();
		TIter next(list) ;
		TKey* key;

		while ( key = (TKey*)next() )
		{
			std::cout << "\tk: " << key->GetName() << ", " << key->GetClassName() << ", TMatrixT<double>" << std::endl;
			if (strcmp(key->GetClassName(), "TMatrixT<double>") == 0 &&
			    (strcmp(key->GetName(), "postfit_cov") == 0 ||
			     strcmp(key->GetName(), "skdetfsi") == 0) )
			{
				std::cout << "\t\tget " << key->GetName() << std::endl;
				matrices.push_back( static_cast<TMatrixT<double>*> (key->ReadObj()) );
			}
		}
	}

	std::cout << "Collected " << matrices.size() << " matrices" << std::endl;
	int cols = 0;
	for (int m = 0; m < matrices.size(); ++m)
		cols += matrices[m]->GetNcols();

	TMatrixD *cov = new TMatrixD(cols, cols);
	TMatrixD *cor = new TMatrixD(cols, cols);
	TMatrixD *ide = new TMatrixD(TMatrixD::kUnit, *cor);

	cols = 0;
	for (int m = 0; m < matrices.size(); ++m)
	{
		int lCols = matrices[m]->GetNcols();
		std::cout << "loading " << lCols << " columns" << std::endl;
		for (int c = 0; c < lCols; ++c)
			for (int r = 0; r < lCols; ++r)
			{
				cov->operator()(cols + r, cols + c) = matrices[m]->operator()(r, c);
				//std::cout << "Out(" << cols + r << ", " << cols + c << ") "
				//	  << mm->operator()(cols + r, cols + c)
				//	  << " = m[" << m << "] at (" << r << ", " << c << ") "
				//	  << matrices[m]->operator()(r, c) << std::endl;
			}

		cols += lCols;
	}

	cols = cov->GetNcols();
	for (int r = 0; r < cols; ++r)
		for (int c = r; c < cols; ++c)
		{
			cor->operator()(r, c) = cov->operator()(r, c) /
				std::sqrt(cov->operator()(r, r) * cov->operator()(c, c) );

			if (r != c && std::abs(cor->operator()(r, c)) - 1 < 1e-6)
			{
				cov->operator()(r, c) *= 0.99999; // = (1-1e-5), to make it non singular
				cor->operator()(r, c) *= 0.99999; // = (1-1e-5), to make it non singular
			}

			cor->operator()(c, r) = cor->operator()(r, c);	//it is symmetric
		}

	TFile out("combinedmatrix.root", "RECREATE");
	cov->Write("covariance");
	cor->Write("correlation");
	ide->Write("identity");
	out.Close();

	return 0;
}
