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
	if (argc < 2)
	{
		std::cerr << "Addmatrix: need at least two parameters\n"
			<< "           output input [input2 ...]\n";
		return 1;
	}

	std::vector<TMatrixT<double>*> matrices;
	matrices.reserve(2*(argc-2));
	std::cout << "Reading " << argc-2 << " files" << std::endl;
	for (int f = 2; f < argc; ++f)
	{
		TFile inf(argv[f], "OPEN");
		if (inf.IsZombie())
			continue;
		std::cout << "opening " << argv[f] << std::endl;

		if (inf.Get("postfit_cov"))
			matrices.push_back( static_cast<TMatrixT<double>*> (inf.Get("postfit_cov")) );
		if (inf.Get("skdetfsi"))
			matrices.push_back( static_cast<TMatrixT<double>*> (inf.Get("skdetfsi")) );


		//TList *list = inf.GetListOfKeys();
		//TIter next(list) ;
		//TKey* key;

		//while ( key = (TKey*)next() )
		//{
		//	std::cout << "\tk: " << key->GetName() << ", " << key->GetClassName() << ", TMatrixT<double>" << std::endl;
		//	if (strcmp(key->GetClassName(), "TMatrixT<double>") == 0 &&
		//	    (strcmp(key->GetName(), "postfit_cov") == 0 ||
		//	     strcmp(key->GetName(), "skdetfsi") == 0) )
		//	{
		//		std::cout << "\t\tget " << key->GetName() << std::endl;
		//		matrices.push_back( static_cast<TMatrixT<double>*> (key->ReadObj()) );
		//	}
		//}
	}

	std::cout << "Collected " << matrices.size() << " matrices" << std::endl;

	if (!matrices.size())
		throw std::invalid_argument("Input files do not contain matrices\n");

	int cols = std::accumulate(matrices.begin(), matrices.end(), 0,
			[](const int &sum, TMatrixT<double> *m) { return sum + m->GetNcols() ;});

	TMatrixD *cov = new TMatrixD(cols, cols);
	TMatrixD *cor = new TMatrixD(cols, cols);
	TMatrixD *ide = new TMatrixD(TMatrixD::kUnit, *cor);

	cols = 0;
	for (const auto &m : matrices) {
		int lCols = m->GetNcols();
		std::cout << "loading " << lCols << " columns" << std::endl;
		for (int c = 0; c < lCols; ++c)
			for (int r = 0; r < lCols; ++r)
			{
				(*cov)(cols + r, cols + c) = (*m)(r, c);
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
			(*cor)(r, c) = (*cov)(r, c) /
				std::sqrt((*cov)(r, r) * (*cov)(c, c) );

			if (r != c && std::abs((*cor)(r, c)) - 1 < 1e-6)
			{
				(*cov)(r, c) *= 0.99999; // = (1-1e-5), to make it non singular
				(*cor)(r, c) *= 0.99999; // = (1-1e-5), to make it non singular
			}

			(*cor)(c, r) = (*cor)(r, c);	//it is symmetric
		}

	TFile out(argv[1], "RECREATE");
	cov->Write("covariance");
	cor->Write("correlation");
	ide->Write("identity");
	out.Close();

	return 0;
}
