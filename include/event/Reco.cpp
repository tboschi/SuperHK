#include "Reco.h"

//ctor
Reco::Reco(std::string cardReco)
{
	CardDealer cd(cardReco);
	std::string fileReco;
	double weight;
	bool gt = cd.Get("reco_path", fileReco);
	if (!cd.Get("scale", weight))
		weight = 1.0;

	//std::cout << "opening " << fileReco << std::endl;
	TFile* inFile = new TFile(fileReco.c_str(), "READ");
	if (inFile->IsZombie())
		std::cerr << "WARNING - Reco: file " << fileReco
			  << " does not exist" << std::endl;

	std::map<std::string, std::string> mSD;
	//look if there is a particular event
	if (cd.Get("Ring_", mSD))
		for (const auto &is : mSD) {
			//std::cout << "histogram " << is.first << " - " << is.second << std::endl;
			// TH2D is made of (xbins + 2) ( ybins + 2)
			// y0 : x0 x1 x2 x3 x4 ...
			// y1 : x0 x1 x2 x3 x4 ...
			// y2 : x0 x1 x2 x3 x4 ...
			// y3 : x0 x1 x2 x3 x4 ...
			// y4 : x0 x1 x2 x3 x4 ...
			// .....
			// as there are overflow and underflow bins
			// correct array starts from 1 to ybins + 1
			if (!inFile->Get(is.second.c_str()))
				continue;

			TH2D *h2 = static_cast<TH2D*>
				(inFile->Get(is.second.c_str()));
		
			int xs = h2->GetNbinsX();// E_true and n_cols
			int ys = h2->GetNbinsY();// E_reco and n_rows
			// catching xs+2 X y2 section
			Eigen::MatrixXd rm = Eigen::Map<Eigen::MatrixXd>
				(h2->GetArray() + xs + 2, xs + 2, ys);
			//std::cout << "raw\n";
			//for (int g = 0; g < 801; ++g)
			//	std::cout << *(h2->GetArray()+g) << "\t";
			//std::cout << std::endl;

			// remove first and last columns
			rm.transposeInPlace();
			rm.leftCols(xs) = rm.middleCols(1, xs);
			rm.conservativeResize(ys, xs);

			rm *= weight;
			_reco[is.first] = rm;

			const double *bx = h2->GetXaxis()->GetXbins()->GetArray();
			const double *by = h2->GetYaxis()->GetXbins()->GetArray();
			_binX[is.first].assign(bx, bx + xs + 1);
			_binY[is.first].assign(by, by + ys + 1);
		}

	inFile->Close();
}

// happy with default copy constructor and destructor

//copy ctor
//Reco::Reco(const Reco & r)
//{
//}

//detor
//Reco::~Reco()
//{
//}


//Get functions

Eigen::MatrixXd Reco::operator()(std::string name)
{
	if (_reco.count(name))
		return _reco[name];
	else {
		std::cerr << "WARNING - Reco: the reconstruction matrix \""
			  << name << "\" does not exist\n";
		return _reco.begin()->second;
	}
}

std::vector<double> Reco::BinsX(std::string name)
{
	if (_binX.count(name))
		return _binX[name];
	else
		return _binX.begin()->second;
}

std::vector<double> Reco::BinsY(std::string name)
{
	if (_binY.count(name))
		return _binY[name];
	else
		return _binY.begin()->second;
}

// rescale all matrices
void Reco::Scale(double X)
{
	for (auto irh = _reco.begin(); irh != _reco.end(); ++irh)
		Scale(irh->first, X);
}

// rescale matrix by name
void Reco::Scale(std::string name, double X)
{
	if (_reco.count(name))
		_reco[name] *= X;
	else
		std::cerr << "WARNING - Reco: the reconstruction matrix \""
			  << name << "\" does not exist\n";
}

// apply vector to matrix to do both scaling and project
// axis is the TH1D direction, so if axis = 'x', then result is in 'y' direction
Eigen::MatrixXd Reco::Apply(std::string name, TH1D* h, char axis)
{
	Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>
		(h->GetArray()+1, h->GetNbinsX(), 1);
	return Apply(name, vec, axis);
}

Eigen::MatrixXd Reco::Apply(std::string name, Eigen::VectorXd &vec, char axis)
{
	switch (axis) {
		case 'x':
			return _reco[name] * vec;
		case 'y':
			return _reco[name].transpose() * vec;
	}
}
